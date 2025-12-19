###############################################################################
# Partial correlations between cognitive scores (MOCA) and pituitary subregional
# quantitative measures across diagnostic groups:
#
# - Computes Spearman partial correlations controlling for Age and Sex
# - Exports results to CSV
# - Generates a block-split square heatmap (ComplexHeatmap)
# - Generates single-page scatterplots for significant associations
#
# NOTE:
# - Users must configure file paths below before execution.
# - No machine-specific absolute paths remain in this script.
# - Intended for open release in GitHub repositories.
###############################################################################

suppressPackageStartupMessages({
  libs <- c(
    "readxl","dplyr","tidyr","tibble","stringr","readr",
    "ppcor",
    "ComplexHeatmap","circlize","grid",
    "ggplot2","ggpubr","ggprism","ggExtra",
    "ragg","glue"
  )
  need <- libs[!libs %in% installed.packages()[,"Package"]]
  if(length(need)) install.packages(need, dependencies = TRUE)
  invisible(lapply(libs, library, character.only = TRUE))
})

###############################################################################
# USER CONFIGURATION
###############################################################################

# Required input Excel file with variables:
# Diagnosis, Age, Sex, MOCA, anterior_Quant_data, posterior_Quant_data,
# stalk_Quant_data, total_Quant_data
infile  <- "path/to/new_final_list1_update.xlsx"   # <-- set manually

# Output root directory (recommended: project outputs folder)
out_dir <- "outputs/partial_corr_results"          # <-- set manually
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Output artifact paths
out_csv <- file.path(out_dir, "partial_corr_MOCA_by_dx.csv")
out_png <- file.path(out_dir, "partial_corr_MOCA_square_split.png")
out_pdf <- file.path(out_dir, "partial_corr_MOCA_square_split.pdf")
scatter_dir <- file.path(out_dir, "scatterplots")

# Analysis parameters
r_min <- -0.25
r_max <-  0.25
alpha_thr <- 0.05   # significance threshold

###############################################################################
# LOAD DATA
###############################################################################
dat <- readxl::read_excel(infile)

need_cols <- c(
  "Diagnosis","MOCA","Age","Sex",
  "anterior_Quant_data","posterior_Quant_data",
  "stalk_Quant_data","total_Quant_data"
)
miss <- setdiff(need_cols, names(dat))
if(length(miss)>0){
  stop("Missing required columns: ", paste(miss, collapse=", "))
}

###############################################################################
# DATA PREP
###############################################################################

dx_levels <- c("MCI","AD","PD","SVD","MS","AQP4Pos_NMOSD")

normalize_dx <- function(x){
  x2 <- stringr::str_trim(as.character(x))
  x2 <- stringr::str_replace_all(x2, fixed("AQP4⁺"), "AQP4+")
  x2
}

coerce_sex01 <- function(x){
  if(is.numeric(x)) return(ifelse(is.na(x), NA_real_, ifelse(x>0,1,0)))
  x2 <- tolower(trimws(as.character(x)))
  ifelse(x2 %in% c("f","female","女"), 1,
         ifelse(x2 %in% c("m","male","男"), 0, NA_real_))
}

dat <- dat %>%
  mutate(
    Diagnosis = factor(normalize_dx(Diagnosis), levels = dx_levels),
    Sex01     = coerce_sex01(Sex),
    Age       = suppressWarnings(as.numeric(Age)),
    MOCA      = suppressWarnings(as.numeric(MOCA)),
    across(ends_with("_Quant_data"), as.numeric)
  ) %>%
  filter(!is.na(Diagnosis))

struct_vars <- c("anterior_Quant_data","posterior_Quant_data",
                 "stalk_Quant_data","total_Quant_data")
struct_lab  <- c("Anterior","Posterior","Stalk","Total")

r_mat <- matrix(NA_real_, nrow = length(struct_vars), ncol = length(dx_levels),
                dimnames=list(struct_lab, dx_levels))
p_mat <- r_mat

###############################################################################
# PARTIAL CORRELATION LOOP
###############################################################################

for(dx in dx_levels){
  dsub <- dat %>% filter(Diagnosis == dx)
  for(k in seq_along(struct_vars)){
    v <- struct_vars[k]
    tmp <- dsub %>% select(MOCA, all_of(v), Age, Sex01) %>% na.omit()
    if(nrow(tmp)>=5){
      pc <- try(ppcor::pcor.test(tmp$MOCA, tmp[[v]],
                                 tmp %>% select(Age,Sex01),
                                 method="spearman"), silent=TRUE)
      if(!inherits(pc,"try-error")){
        r_mat[k,dx] <- unname(pc$estimate)
        p_mat[k,dx] <- unname(pc$p.value)
      }
    }
  }
}

###############################################################################
# EXPORT TIDY RESULTS
###############################################################################

res_tbl <- as.data.frame(r_mat) %>%
  tibble::rownames_to_column("Structure") %>%
  pivot_longer(-Structure, names_to="Diagnosis", values_to="partial_r")

p_tbl <- as.data.frame(p_mat) %>%
  tibble::rownames_to_column("Structure") %>%
  pivot_longer(-Structure, names_to="Diagnosis", values_to="p_value")

out_tbl <- left_join(res_tbl,p_tbl,by=c("Structure","Diagnosis"))
readr::write_csv(out_tbl, out_csv)

###############################################################################
# HEATMAP (no on-screen preview)
###############################################################################

pal <- colorRampPalette(c("#1F4E79","#6AAED6","white","#E57373","#B22222"))(101)
col_fun <- circlize::colorRamp2(seq(r_min,r_max,length.out=length(pal)), pal)

dx_display <- c("MCI","AD","PD","SVD","MS","NMOSD")
colnames(r_mat)[colnames(r_mat)=="AQP4Pos_NMOSD"] <- "NMOSD"
colnames(p_mat)[colnames(p_mat)=="AQP4Pos_NMOSD"] <- "NMOSD"

dx_colors <- c("MCI"="#C62E2E","AD"="#F05A24","PD"="#E2A019",
               "SVD"="#58A6E8","MS"="#2F8F46","NMOSD"="#D895C0")

top_anno <- HeatmapAnnotation(
  Dx = factor(colnames(r_mat), levels=dx_display),
  col = list(Dx = dx_colors),
  show_legend=FALSE,
  height=unit(6,"mm")
)

sig_symbol <- function(p){
  if (is.na(p)) "" else if (p<0.001) "***"
  else if (p<0.01) "**" else if (p<0.05) "*" else ""
}

col_split <- factor(colnames(r_mat), levels=dx_display)

ht <- Heatmap(
  r_mat,
  col = col_fun,
  name="Partial_r",
  top_annotation=top_anno,
  cluster_rows=FALSE, cluster_columns=FALSE,
  column_split = col_split,
  gap = unit(4,"mm"),
  row_names_side="right",
  show_column_names=FALSE,
  na_col = "#F7F7F7",
  cell_fun = function(j,i,x,y,w,h,fill){
    lab <- sig_symbol(p_mat[i,j])
    if(nzchar(lab)) grid.text(lab,x,y,gp=gpar(fontsize=12,fontface="bold"))
  }
)

ragg::agg_png(out_png,width=1700,height=1350,units="px",res=300)
draw(ht, merge_legends = TRUE)
dev.off()

pdf(out_pdf,width=6.6,height=5.2,useDingbats=FALSE)
draw(ht, merge_legends = TRUE)
dev.off()

###############################################################################
# SIGNIFICANT SCATTERPLOTS
###############################################################################

sig_idx <- which(p_mat <= (alpha_thr+1e-12), arr.ind=TRUE)

if(nrow(sig_idx)>0){
  
  if(!dir.exists(scatter_dir)) dir.create(scatter_dir,recursive=TRUE)
  
  struct_to_var <- setNames(struct_vars, struct_lab)
  display_to_data <- c("NMOSD"="AQP4Pos_NMOSD",
                       setNames(dx_levels, dx_levels))
  
  get_rank_resid <- function(df,y_var,x_var="MOCA"){
    d2 <- df[,c(x_var,y_var,"Age","Sex01")] |> na.omit()
    if(nrow(d2)<5) return(NULL)
    rx <- residuals(lm(rank(d2[[x_var]])~Age+Sex01,data=d2))
    ry <- residuals(lm(rank(d2[[y_var]])~Age+Sex01,data=d2))
    data.frame(x_resid=rx,y_resid=ry)
  }
  
  idx_rows <- list()
  
  for (k in seq_len(nrow(sig_idx))){
    i <- sig_idx[k,1]; j <- sig_idx[k,2]
    struct_disp <- rownames(p_mat)[i]
    dx_disp     <- colnames(p_mat)[j]
    dx_data     <- display_to_data[[dx_disp]]
    y_var <- struct_to_var[[struct_disp]]
    
    dsub <- dat %>%
      filter(Diagnosis==dx_data) %>%
      select(MOCA, all_of(struct_vars), Age, Sex01)
    
    df_plot <- get_rank_resid(dsub,y_var=y_var)
    if(is.null(df_plot)) next
    
    ct <- suppressWarnings(cor.test(df_plot$x_resid,df_plot$y_resid,
                                    method="pearson"))
    
    base <- glue("scatter_{dx_disp}_{struct_disp}")
    png_file <- file.path(scatter_dir,paste0(base,".png"))
    pdf_file <- file.path(scatter_dir,paste0(base,".pdf"))
    
    p_sc <- ggplot(df_plot,aes(x=x_resid,y=y_resid)) +
      geom_point() +
      geom_smooth(method="lm",se=FALSE) +
      ggpubr::stat_cor(method="pearson") +
      ggprism::theme_prism(border=TRUE) +
      labs(title=glue("{dx_disp}: {struct_disp} vs. MOCA (partial Spearman via ranks)"))
    
    pm <- ggExtra::ggMarginal(p_sc,type="histogram",margins="both")
    
    ragg::agg_png(png_file,width=1200,height=900,units="px",res=150)
    grid.draw(pm); dev.off()
    
    pdf(pdf_file,width=8,height=6,useDingbats=FALSE)
    grid.draw(pm); dev.off()
    
    idx_rows[[length(idx_rows)+1]] <- data.frame(
      Diagnosis_display=dx_disp,
      Diagnosis_data=dx_data,
      Structure=struct_disp,
      r=unname(ct$estimate),
      p=unname(ct$p.value),
      n=nrow(df_plot),
      png=png_file,
      pdf=pdf_file
    )
  }
  
  if(length(idx_rows)){
    write_csv(bind_rows(idx_rows), file.path(scatter_dir,"scatter_index.csv"))
  }
}

message("Partial correlation pipeline completed.")