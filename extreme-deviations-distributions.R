
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(fmsb)
})

# -------- 数据读取与处理 -----------------------------------------------------
csv_path <- 'E:/Lifespan_freesurfer_results/Github/Test_results/V7_DK/extreme_deviation_positive.csv'
dat <- read_csv(csv_path, show_col_types = FALSE)

need <- c("category","anterior_pituitary","posterior_pituitary",
          "pituitary_stalk","total_pituitary")
stopifnot(all(need %in% names(dat)))

# 固定疾病组顺序
dx_order <- c("MCI","AD","PD","SVD","MS","NMOSD")
dat$category <- factor(dat$category, levels = dx_order)

# 转置：行 = 结构，列 = 疾病组
long <- dat %>% pivot_longer(-category, names_to = "structure", values_to = "value")
wide <- long %>% pivot_wider(names_from = category, values_from = value)
mat  <- as.data.frame(wide); rownames(mat) <- mat$structure; mat <- mat[, -1]
mat  <- mat[, dx_order, drop = FALSE]

# ★ 绘制顺序：Stalk -> Anterior -> Posterior -> Total（最后一个在最上层）
draw_order <- c("pituitary_stalk", "anterior_pituitary",
                "posterior_pituitary", "total_pituitary")
mat <- mat[draw_order, , drop = FALSE]

# fmsb 需要 max/min 两行
vmax  <- max(mat, na.rm = TRUE)
y_top <- max(0.25, ceiling(vmax/0.05)*0.05)
plot_df <- rbind(rep(y_top, ncol(mat)), rep(0, ncol(mat)), mat)
rownames(plot_df) <- c("max","min", rownames(mat))

# ★ 颜色顺序与行顺序一一对应：Stalk(红) -> Anterior(天蓝) -> Posterior(深蓝) -> Total(黄)
pcol  <- c("#CA0E12", "#2AA7DE", "#25377F", "#F6BD21")
pfcol <- sapply(pcol, function(z){
  rgb_val <- grDevices::col2rgb(z)/255
  do.call(rgb, as.list(c(rgb_val, 0.35)))
})

# -------- 导出 PDF -----------------------------------------------------------
pdf("/Users/bailu/Desktop/radar_plot.pdf", width = 6, height = 6)

op <- par(mar = c(2.2, 2, 2, 2)); par(xpd = NA)
fmsb::radarchart(
  plot_df,
  axistype = 1, seg = 5, caxislabels = rep("", 5),
  cglcol = "grey70", cglty = 3, cglwd = 0.8,
  vlcex = 1.05, pcol = pcol, pfcol = pfcol, plwd = 2, plty = 1
)
legend("bottom", inset = c(0, -0.06),
       legend = rownames(plot_df)[-(1:2)],
       col = pcol, pch = 16, pt.cex = 1.6, cex = 0.95,
       bty = "n", ncol = 4, xpd = NA, y.intersp = 1.2, x.intersp = 0.8)
par(op); dev.off()

# -------- 导出 PNG -----------------------------------------------------------
png("csv_path <- 'E:/Lifespan_freesurfer_results/Github/Test_results/V7_DK/radar_plot.png", width = 2000, height = 2000, res = 300)

op <- par(mar = c(2.2, 2, 2, 2)); par(xpd = NA)
fmsb::radarchart(
  plot_df,
  axistype = 1, seg = 5, caxislabels = rep("", 5),
  cglcol = "grey70", cglty = 3, cglwd = 0.8,
  vlcex = 1.05, pcol = pcol, pfcol = pfcol, plwd = 2, plty = 1
)
legend("bottom", inset = c(0, -0.06),
       legend = rownames(plot_df)[-(1:2)],
       col = pcol, pch = 16, pt.cex = 1.6, cex = 0.95,
       bty = "n", ncol = 4, xpd = NA, y.intersp = 1.2, x.intersp = 0.8)
par(op); dev.off()