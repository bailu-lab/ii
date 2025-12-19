###############################################################################
# Survival analyses for pituitary metrics:
# Cox regression, Kaplan–Meier curves, and forest plots
#
# Part A: Conversion cohorts (HC→MCI, MCI→AD)
# Part B: MS & NMOSD cohorts (three pituitary structures, raw volumes)
#
# This script is intended for public release on GitHub.
# Users are expected to supply data files and configure IO paths below.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(broom)
  library(readr)
  library(grid)
  library(gtable)
})

## ---------------------------------------------------------------------------
## USER CONFIGURATION
## ---------------------------------------------------------------------------
# provide paths on your machine before running
infile_conv  <- "path/to/longitudinal_dataset.xlsx"
infile_ms    <- "path/to/disease_followup_dataset.xlsx"

out_dir_conv <- "outputs/cox_km_conversion"
out_dir_ms   <- "outputs/cox_km_MS_NMOSD"

dir.create(out_dir_conv, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_ms,  recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## common utilities
## ---------------------------------------------------------------------------

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

km_palette <- c("#8CA5EA", "#E3B13C")   # low / high risk

plot_km_unified <- function(fit_km, data, strata_labels, palette,
                            time_scale   = 1,
                            xlab         = "Follow-up (years)",
                            ylab         = "Survival probability",
                            title        = NULL,
                            subtitle     = NULL,
                            vline_median = NA_real_,
                            x_start0     = 0,
                            p_label      = NULL) {
  
  ss <- survminer::surv_summary(fit_km, data = data)
  ss$time_scaled <- ss$time * time_scale
  
  if ("strata" %in% names(ss)) {
    ss$strata <- factor(ss$strata)
    levels(ss$strata) <- strata_labels
  }
  
  ## insert t=0
  if ("strata" %in% names(ss)) {
    levs <- levels(ss$strata)
    add0 <- data.frame(
      time        = 0, time_scaled = 0, surv = 1,
      lower = 1, upper = 1,
      strata = factor(levs, levels = levs)
    )
    ss <- rbind(ss[, intersect(names(ss), names(add0))], add0)
    ss <- ss[order(ss$strata, ss$time_scaled), ]
  } else {
    ss <- rbind(
      data.frame(time=0,time_scaled=0,surv=1,lower=1,upper=1),
      ss
    )
  }
  
  censor_df <- NULL
  if ("n.censor" %in% names(ss)) {
    idx <- which(!is.na(ss$n.censor) & ss$n.censor > 0)
    if (length(idx)) censor_df <- ss[idx, ]
  }
  
  xmax <- max(ss$time_scaled, na.rm = TRUE)
  
  g <- ggplot(ss, aes(x = time_scaled, y = surv, color = strata)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = strata),
                alpha = 0.14, color = NA, show.legend = FALSE) +
    geom_step(linewidth = 1.4) +
    scale_color_manual(values = palette, labels = strata_labels, name = NULL) +
    scale_fill_manual(values = palette, guide = "none") +
    scale_y_continuous(limits = c(0, 1.02)) +
    scale_x_continuous(limits = c(x_start0, xmax)) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      legend.position = "top"
    )
  
  if (!is.null(censor_df)) {
    g <- g + geom_point(data = censor_df,
                        aes(x = time_scaled, y = surv, color = strata),
                        shape = 3, size = 2.4)
  }
  
  if (is.finite(vline_median)) {
    g <- g + geom_vline(xintercept = vline_median,
                        linetype = "dashed", linewidth = 0.6)
  }
  
  if (!is.null(p_label)) {
    g <- g + annotate(
      "text",
      x = x_start0 + 0.05 * (xmax - x_start0),
      y = 0.18,
      label = p_label,
      hjust = 0, size = 5
    )
  }
  
  return(g)
}

## ===========================================================================
## PART A: conversion cohorts (HC→MCI, MCI→AD)
## ===========================================================================

run_cox_conversion <- function(infile_path, out_dir){
  
  dat_raw <- readxl::read_xlsx(infile_path) |> as.data.frame()
  names(dat_raw) <- trimws(names(dat_raw))
  
  req_cols <- c("baseline_for_follow_up","Diagnosis",
                "follow_up_years",
                "anterior_Quant_data",
                "posterior_Quant_data",
                "stalk_Quant_data")
  if (!all(req_cols %in% names(dat_raw))) {
    stop("Missing required columns for conversion cohorts.")
  }
  
  # clean numeric
  conv0 <- dat_raw |>
    mutate(
      SubjectID = baseline_for_follow_up,
      Diagnosis = as.character(Diagnosis),
      follow_up_years = to_num(follow_up_years),
      anterior_Quant_data  = to_num(anterior_Quant_data),
      posterior_Quant_data = to_num(posterior_Quant_data),
      stalk_Quant_data     = to_num(stalk_Quant_data)
    ) |>
    filter(!is.na(SubjectID),
           Diagnosis %in% c("HC","MCI","AD"),
           is.finite(follow_up_years)) |>
    arrange(SubjectID, follow_up_years)
  
  base_df <- conv0 |>
    group_by(SubjectID) |>
    slice_min(order_by = follow_up_years, n = 1, with_ties = FALSE)
  
  dat <- conv0 |>
    left_join(base_df |> transmute(SubjectID,
                                   base_time = follow_up_years,
                                   base_diag = Diagnosis),
              by = "SubjectID")
  
  build_cohort <- function(data, baseline_diag, target_diag){
    ids <- data |> filter(Diagnosis == baseline_diag,
                          follow_up_years == base_time) |> pull(SubjectID)
    df <- data |> filter(SubjectID %in% ids)
    
    conv <- df |>
      group_by(SubjectID) |>
      summarise(
        t0   = min(follow_up_years),
        tmax = max(follow_up_years),
        tevt = suppressWarnings(min(follow_up_years[Diagnosis == target_diag])),
        .groups = "drop"
      ) |>
      mutate(event = ifelse(is.infinite(tevt), 0, 1),
             time  = ifelse(event==1, tevt - t0, tmax - t0))
    
    base_cov <- df |>
      filter(follow_up_years == base_time) |>
      distinct(SubjectID,.keep_all = TRUE) |>
      transmute(
        SubjectID,
        anterior_Quant_data,
        posterior_Quant_data,
        stalk_Quant_data
      )
    
    res <- left_join(conv, base_cov, by = "SubjectID") |>
      filter(time > 0)
    attr(res,"desc") <- paste0(baseline_diag,"→",target_diag)
    return(res)
  }
  
  coh1 <- build_cohort(dat,"HC","MCI")
  coh2 <- build_cohort(dat,"MCI","AD")
  cohorts <- list(coh1, coh2)
  
  for (coh in cohorts) {
    nm <- attr(coh,"desc")
    if (nrow(coh) < 10) next
    
    fml <- as.formula("Surv(time,event) ~ anterior_Quant_data + posterior_Quant_data + stalk_Quant_data")
    fit <- coxph(fml, data = coh, x = TRUE, y = TRUE)
    
    coh$lp <- as.numeric(predict(fit,type="lp"))
    thr <- median(coh$lp,na.rm=TRUE)
    coh$risk_grp <- factor(ifelse(coh$lp > thr,"High risk","Low risk"),
                           levels=c("Low risk","High risk"))
    
    fit_km <- survfit(Surv(time,event) ~ risk_grp, data = coh)
    sd_obj <- survdiff(Surv(time,event) ~ risk_grp, data = coh)
    df_chi <- max(1,length(sd_obj$n)-1)
    p_lr <- pchisq(sd_obj$chisq,df=df_chi,lower.tail=FALSE)
    p_txt <- if (p_lr<0.001) "p < 0.001" else sprintf("p = %.3f",p_lr)
    
    tab_km <- summary(fit_km)$table
    vline_med <- if ("median" %in% colnames(tab_km))
      unname(tab_km[2,"median"]) else NA_real_
    
    g_km <- plot_km_unified(
      fit_km,data=coh,
      strata_labels=levels(coh$risk_grp),
      palette=km_palette,
      title=paste0("KM: ",nm),
      subtitle="Risk split by Cox linear predictor (median)",
      vline_median=vline_med
    )
    
    base <- file.path(out_dir,paste0("KM_",gsub("→","to_",nm)))
    ggsave(paste0(base,".png"),g_km,width=8,height=6,dpi=300,bg="white")
    ggsave(paste0(base,".pdf"),g_km,width=8,height=6,bg="white")
    
    tab <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
      filter(term %in% c("anterior_Quant_data","posterior_Quant_data","stalk_Quant_data"))
    
    write_csv(tab |> mutate(Cohort = nm),
              file.path(out_dir,paste0("Cox_terms_",gsub("→","to_",nm),".csv")))
  }
}

## ===========================================================================
## PART B: MS / NMOSD
## ===========================================================================

run_cox_ms <- function(infile_path,out_dir,time_var="Followup_time_month"){
  
  raw <- readxl::read_xlsx(infile_path) |> as.data.frame()
  names(raw) <- trimws(names(raw))
  
  struct_cols <- c("anterior_Quant_data","posterior_Quant_data","stalk_Quant_data")
  if (!all(struct_cols %in% names(raw))) stop("missing pituitary columns")
  
  df <- raw |>
    transmute(
      Diagnosis = as.character(Diagnosis),
      t = to_num(.data[[time_var]]),
      pit_anterior  = to_num(anterior_Quant_data),
      pit_posterior = to_num(posterior_Quant_data),
      pit_stalk     = to_num(stalk_Quant_data),
      EDSS_prog     = to_num(EDSS_Progress),
      SPMS_conv     = to_num(SPMS_conversion),
      relapse       = to_num(Followup_relapse)
    )
  
  outcomes <- c("EDSS_prog","SPMS_conv","relapse")
  diseases <- intersect(c("MS","AQP4Pos_NMOSD"), unique(df$Diagnosis))
  
  for (diag in diseases) {
    sub <- df |> filter(Diagnosis==diag)
    
    for (ev in outcomes){
      if (!ev %in% names(sub)) next
      
      sub1 <- sub |> filter(!is.na(t), t>0)
      if (nrow(sub1)<15) next
      
      evt <- ifelse(sub1[[ev]]>0,1,0)
      if (sum(evt)<5) next
      
      fml <- as.formula(sprintf("Surv(t, evt) ~ pit_anterior+pit_posterior+pit_stalk"))
      fit <- coxph(fml,data=sub1,x=TRUE,y=TRUE)
      
      sub1$lp <- as.numeric(predict(fit,type="lp"))
      thr <- median(sub1$lp,na.rm=TRUE)
      sub1$risk_grp <- factor(ifelse(sub1$lp>thr,"High risk","Low risk"),
                              levels=c("Low risk","High risk"))
      
      fit_km <- survfit(Surv(t,evt) ~ risk_grp,data=sub1)
      sd_obj <- survdiff(Surv(t,evt) ~ risk_grp,data=sub1)
      df_chi <- max(1,length(sd_obj$n)-1)
      p_lr <- pchisq(sd_obj$chisq,df=df_chi,lower.tail=FALSE)
      p_txt <- if (p_lr<0.001) "p < 0.001" else sprintf("p = %.3f",p_lr)
      
      tab_km <- summary(fit_km)$table
      vline_med <- if ("median" %in% colnames(tab_km))
        unname(tab_km[2,"median"])/12 else NA_real_
      
      g_km <- plot_km_unified(
        fit_km,data=sub1,
        strata_labels=levels(sub1$risk_grp),
        palette=km_palette,
        time_scale=1/12,
        title=paste0("KM survival: ",ev),
        subtitle=paste0("Diagnosis: ",diag),
        vline_median=vline_med,
        p_label=p_txt
      )
      
      base <- file.path(out_dir,paste0("KM_",ev,"_",diag))
      ggsave(paste0(base,".png"),g_km,width=8,height=6,dpi=300,bg="white")
      ggsave(paste0(base,".pdf"),g_km,width=8,height=6,bg="white")
      
      tab <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
        filter(term %in% c("pit_anterior","pit_posterior","pit_stalk"))
      
      write_csv(tab |> mutate(Diagnosis=diag, Outcome=ev),
                file.path(out_dir,paste0("Cox_terms_",ev,"_",diag,".csv")))
    }
  }
}

## ===========================================================================
## EXECUTION
## ===========================================================================

run_cox_conversion(infile_conv,out_dir_conv)
run_cox_ms(infile_ms,out_dir_ms)