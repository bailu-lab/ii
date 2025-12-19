#!/usr/bin/env Rscript
##
## ROC_all_groups_logistic_CV10.R
##
## Purpose:
##   Perform 10-fold cross-validated logistic regression classification
##   for multiple disease groups vs. healthy controls (HC) using
##   pituitary-derived features, and:
##     1) compute pooled ROC curves with DeLong 95% CIs
##     2) export per-fold AUCs and cross-validated confusion matrices
##
## Requirements on input data:
##   - Excel file with at least the following columns:
##       Diagnosis    : factor/character, must contain "HC" and disease labels
##       Age          : numeric
##       Sex          : factor/character (will be coerced to factor)
##       anterior_Quant_data
##       posterior_Quant_data
##       stalk_Quant_data
##       total_Quant_data
##
## Usage (typical):
##   - Place your Excel file under ./data/
##   - Set `infile` below to the relative path, then run:
##       source("scripts/ROC_all_groups_logistic_CV10.R")
##

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(MatchIt)
  library(pROC)
  library(ggplot2)
  library(RColorBrewer)
  library(readr)
})

## -------- 0) User parameters --------

# Input file: relative path within the project
infile  <- "data/real_data.xlsx"

# Output directory (will be created if not existing)
out_dir <- "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Output figure paths
out_png <- file.path(out_dir, "ROC_all_groups_logistic_CV10.png")
out_pdf <- file.path(out_dir, "ROC_all_groups_logistic_CV10.pdf")

# Cross-validation parameters
K_folds <- 10
cv_seed <- 2025

# Disease groups to be compared with HC
# Make sure these labels exist in the Diagnosis column of your data
disease_levels <- c("MCI", "AD", "PD", "SVD", "MS", "AQP4Pos_NMOSD")

# Pretty names for plotting
pretty_name <- function(dx) ifelse(dx == "AQP4Pos_NMOSD", "NMOSD", dx)

## -------- 1) Read data & basic cleaning --------

dat0 <- readxl::read_xlsx(infile)

# Standardize column names (remove spaces, unify style)
names(dat0) <- names(dat0) |>
  stringr::str_replace_all("\\s+", "_") |>
  stringr::str_trim()

req_cols <- c(
  "Diagnosis", "Age", "Sex",
  "anterior_Quant_data", "posterior_Quant_data",
  "stalk_Quant_data", "total_Quant_data"
)
stopifnot(all(req_cols %in% names(dat0)))

dat <- dat0 %>%
  mutate(
    Diagnosis = as.character(Diagnosis),
    Sex       = as.factor(Sex)
  ) %>%
  filter(!is.na(Diagnosis), !is.na(Age), !is.na(Sex)) %>%
  drop_na(any_of(c(
    "anterior_Quant_data", "posterior_Quant_data",
    "stalk_Quant_data", "total_Quant_data"
  )))

## ========================================================
## A) Single matched analysis on the full dataset (reference only)
##    - Not used for plotting / confusion matrices
## ========================================================

fit_one_group <- function(dx, df) {
  cases <- df %>% filter(Diagnosis == dx)
  ctrls <- df %>% filter(Diagnosis == "HC")
  if (nrow(cases) == 0 || nrow(ctrls) == 0) {
    warning(sprintf("Group %s cannot be matched: zero cases or controls.", dx))
    return(NULL)
  }
  tmp <- bind_rows(
    cases %>% mutate(group = 1L),
    ctrls %>% mutate(group = 0L)
  )
  tmp$group <- as.integer(tmp$group)
  
  m  <- MatchIt::matchit(
    group ~ Age + Sex,
    data   = tmp,
    method = "nearest",
    ratio  = 1,
    replace = FALSE
  )
  md <- MatchIt::match.data(m)
  
  fit <- glm(
    group ~ anterior_Quant_data + posterior_Quant_data +
      stalk_Quant_data + total_Quant_data,
    data   = md,
    family = binomial()
  )
  
  pr      <- predict(fit, type = "response")
  roc_obj <- pROC::roc(md$group, pr, quiet = TRUE, ci = TRUE, ci.method = "delong")
  ci_auc  <- pROC::ci.auc(roc_obj)
  
  list(
    dx      = dx,
    pretty  = pretty_name(dx),
    n_case  = sum(md$group == 1),
    n_ctrl  = sum(md$group == 0),
    roc     = roc_obj,
    auc_ci  = ci_auc,
    fit     = fit
  )
}

results_ref <- lapply(disease_levels, fit_one_group, df = dat)
results_ref <- Filter(Negate(is.null), results_ref)

cat("===== Reference (single matched dataset, not used for CV plots) =====\n")
if (length(results_ref) > 0) {
  for (x in results_ref) {
    ci <- as.numeric(x$auc_ci)
    cat(sprintf(
      "%-7s  n_case=%-4d  n_ctrl=%-4d   AUC=%.3f  (%.3f–%.3f)\n",
      x$pretty, x$n_case, x$n_ctrl, ci[2], ci[1], ci[3]
    ))
  }
}
cat("=====================================================================\n")

## ========================================================
## B) 10-fold CV: pooled ROC per disease group
## ========================================================

make_stratified_folds <- function(labels, K = 10, seed = 1) {
  set.seed(seed)
  idx1 <- which(labels == 1L)
  idx0 <- which(labels == 0L)
  f1   <- sample(rep(1:K, length.out = length(idx1)))
  f0   <- sample(rep(1:K, length.out = length(idx0)))
  folds <- vector("list", K)
  for (k in seq_len(K)) {
    val_idx <- c(idx1[f1 == k], idx0[f0 == k])
    folds[[k]] <- val_idx
  }
  folds
}

cv10_one_group <- function(dx, df, K = 10, seed = 1) {
  cases <- df %>% filter(Diagnosis == dx)
  ctrls <- df %>% filter(Diagnosis == "HC")
  if (nrow(cases) == 0 || nrow(ctrls) == 0) return(NULL)
  
  dat_all <- bind_rows(
    cases %>% mutate(group = 1L),
    ctrls %>% mutate(group = 0L)
  ) %>% mutate(group = as.integer(group))
  
  folds <- make_stratified_folds(dat_all$group, K = K, seed = seed)
  
  preds <- numeric(0)
  trues <- integer(0)
  auc_per_fold <- rep(NA_real_, K)
  
  for (k in seq_len(K)) {
    val_idx <- folds[[k]]
    train   <- dat_all[-val_idx, ]
    valid   <- dat_all[ val_idx, ]
    
    # Matching in training folds
    m_tr  <- MatchIt::matchit(
      group ~ Age + Sex, data = train,
      method = "nearest", ratio = 1, replace = FALSE
    )
    tr_md <- tryCatch(MatchIt::match.data(m_tr), error = function(e) NULL)
    if (is.null(tr_md) || nrow(tr_md) < 5) next
    
    fit <- glm(
      group ~ anterior_Quant_data + posterior_Quant_data +
        stalk_Quant_data + total_Quant_data,
      data   = tr_md,
      family = binomial()
    )
    
    # Matching in validation folds (within fold)
    m_va  <- MatchIt::matchit(
      group ~ Age + Sex, data = valid,
      method = "nearest", ratio = 1, replace = FALSE
    )
    va_md <- tryCatch(MatchIt::match.data(m_va), error = function(e) NULL)
    if (is.null(va_md) || nrow(va_md) < 5) next
    
    pr    <- predict(fit, newdata = va_md, type = "response")
    roc_k <- pROC::roc(va_md$group, pr, quiet = TRUE)
    auc_per_fold[k] <- as.numeric(pROC::auc(roc_k))
    
    preds <- c(preds, pr)
    trues <- c(trues, va_md$group)
  }
  
  # Pooled ROC curve across all validation folds
  roc_all <- pROC::roc(trues, preds, quiet = TRUE, ci = TRUE, ci.method = "delong")
  
  list(
    dx        = dx,
    pretty    = pretty_name(dx),
    auc_per_fold = auc_per_fold,
    auc_mean  = mean(auc_per_fold, na.rm = TRUE),
    auc_sd    = sd(auc_per_fold,   na.rm = TRUE),
    roc_all   = roc_all,
    preds_all = preds,
    trues_all = trues
  )
}

cv_results <- lapply(disease_levels, cv10_one_group,
                     df = dat, K = K_folds, seed = cv_seed)
cv_results <- Filter(Negate(is.null), cv_results)
stopifnot(length(cv_results) > 0)

## Export CV summary and per-fold AUCs

cv_summary <- lapply(cv_results, function(x) {
  ci <- as.numeric(pROC::ci.auc(x$roc_all))
  data.frame(
    Group    = x$pretty,
    K        = K_folds,
    AUC_mean = round(x$auc_mean, 4),
    AUC_sd   = round(x$auc_sd,   4),
    AUC_ci_l = round(ci[1], 4),
    AUC_ci_m = round(ci[2], 4),
    AUC_ci_u = round(ci[3], 4),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

cv_folds <- lapply(cv_results, function(x) {
  data.frame(
    Group = x$pretty,
    Fold  = seq_along(x$auc_per_fold),
    AUC   = x$auc_per_fold,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()

readr::write_csv(cv_summary, file.path(out_dir, "ROC_CV10_summary.csv"))
readr::write_csv(cv_folds,   file.path(out_dir, "ROC_CV10_folds.csv"))

message("10-fold cross-validated ROC summary:")
print(cv_summary)

## ========================================================
## C) Plot pooled ROC curves (10-fold CV)
## ========================================================

roc_list_cv <- setNames(lapply(cv_results, `[[`, "roc_all"),
                        sapply(cv_results, function(x) x$pretty))

auc_labels_cv <- sapply(cv_results, function(x) {
  ci <- as.numeric(pROC::ci.auc(x$roc_all))
  sprintf("%s  AUC=%.3f [%.3f–%.3f]",
          x$pretty, ci[2], ci[1], ci[3])
})
names(auc_labels_cv) <- sapply(cv_results, `[[`, "pretty")

task_levels <- c("MCI","AD","PD","SVD","MS","NMOSD")
task_levels <- task_levels[task_levels %in% names(roc_list_cv)]

pal <- RColorBrewer::brewer.pal(length(task_levels), "Set1")
names(pal) <- task_levels
label_vec <- auc_labels_cv[task_levels]; names(label_vec) <- task_levels

g <- pROC::ggroc(roc_list_cv[task_levels], legacy.axes = TRUE, size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              linewidth = 1.0, color = "grey55") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_classic(base_size = 16) +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    axis.line         = element_blank(),
    axis.title        = element_text(face = "bold"),
    axis.text         = element_text(color = "black"),
    legend.position   = c(0.30, 0.03),
    legend.justification = c(0, 0),
    legend.title      = element_blank(),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text       = element_text(size = 12, margin = margin(b = 2))
  ) +
  scale_color_manual(values = pal, labels = label_vec) +
  guides(color = guide_legend(
    override.aes = list(linewidth = 1.6),
    keyheight = unit(10, "pt"), keywidth = unit(18, "pt"),
    ncol = 1, byrow = TRUE
  ))

ggsave(out_png, g, width = 5.2, height = 5.2, dpi = 300, bg = "white")
ggsave(out_pdf, g, width = 5.2, height = 5.2)
message("ROC figure saved to:")
message(" - ", out_png)
message(" - ", out_pdf)

## ========================================================
## D) Confusion matrices based on pooled CV predictions
## ========================================================

cm_from_cvresult <- function(x) {
  dx_name  <- x$pretty
  trues    <- as.integer(x$trues_all)   # 0 = HC, 1 = disease
  preds    <- as.numeric(x$preds_all)
  
  thr <- as.numeric(pROC::coords(
    x$roc_all, x = "best",
    best.method = "youden",
    ret = "threshold", transpose = TRUE
  )[1])
  
  pred_label <- ifelse(preds >= thr, 1L, 0L)
  
  TP <- sum(pred_label == 1L & trues == 1L, na.rm = TRUE)
  TN <- sum(pred_label == 0L & trues == 0L, na.rm = TRUE)
  FP <- sum(pred_label == 1L & trues == 0L, na.rm = TRUE)
  FN <- sum(pred_label == 0L & trues == 1L, na.rm = TRUE)
  
  ACC <- (TP + TN) / (TP + TN + FP + FN)
  SEN <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  SPE <- ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_)
  N   <- TP + TN + FP + FN
  
  data.frame(
    True  = factor(c("HC", "HC", dx_name, dx_name), levels = c("HC", dx_name)),
    Pred  = factor(c("HC", dx_name, "HC", dx_name), levels = c("HC", dx_name)),
    Freq  = as.numeric(c(TN, FP, FN, TP)),
    Group = dx_name,
    thr   = as.numeric(thr),
    ACC   = as.numeric(ACC),
    SEN   = as.numeric(SEN),
    SPE   = as.numeric(SPE),
    N     = as.numeric(N),
    check.names = FALSE
  )
}

cm_cv_df <- dplyr::bind_rows(lapply(cv_results, cm_from_cvresult))
stopifnot(nrow(cm_cv_df) > 0)

panel_titles <- cm_cv_df %>%
  dplyr::distinct(Group, thr, N) %>%
  dplyr::mutate(
    Panel = sprintf(
      "%s vs HC  (CV10, thr=%.3f, n=%d, matched-in-folds)",
      Group, thr, N
    )
  )
cm_cv_df <- dplyr::left_join(cm_cv_df, panel_titles, by = c("Group", "thr", "N"))

reds <- c("#fee3ce", "#eabaa1", "#dc917b", "#d16d5b", "#c44438", "#b7282e")

g_cm <- ggplot(cm_cv_df, aes(x = Pred, y = True, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Freq), size = 6.5, fontface = "bold", color = "black") +
  scale_fill_gradientn(colours = reds, guide = "none") +
  labs(x = "Predicted label", y = "True label") +
  facet_wrap(~ Panel, ncol = 3, scales = "free") +
  theme_classic(base_size = 17) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_blank(),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "black"),
    axis.line        = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.spacing    = unit(0.9, "lines"),
    legend.position  = "none"
  )

cm_png <- file.path(out_dir, "CM_disease_vs_HC_matched_CV10.png")
cm_pdf <- file.path(out_dir, "CM_disease_vs_HC_matched_CV10.pdf")
ggsave(cm_png, g_cm, width = 8, height = 5, dpi = 300, bg = "white")
ggsave(cm_pdf, g_cm, width = 8, height = 5)

cm_csv <- file.path(out_dir, "CM_disease_vs_HC_matched_CV10.csv")
readr::write_csv(cm_cv_df, cm_csv)

message("Confusion matrices (10-fold CV pooled predictions) saved to:")
message(" - ", cm_png)
message(" - ", cm_pdf)
message(" - ", cm_csv)