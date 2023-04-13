# This script is designed for final plotting with heatmap
# Weihua Guo
# 08/24/2019
rm(list = ls())


suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))

download_dir <- "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir <- paste(download_dir, "Fig2", sep = "/")
out_dir <- paste(download_dir, "Fig2", sep = "/")


low_corr_file = "_corr_coef_low_pearson_immune_only_xself_.csv"
high_corr_file = "_corr_coef_high_pearson_immune_only_xself_.csv"

cat("Reading correlation coefficients...\n")
low_corr = read.csv(paste(work_dir, low_corr_file, sep = ""), row.names = 1)
high_corr = read.csv(paste(work_dir, high_corr_file, sep = ""), row.names = 1)

st = Sys.time()
cat("Ploting heatmaps (orderd by hot cohort)... \n")
col_fun = colorRamp2(c(-1,0,1), c("darkorange", "cornsilk", "forestgreen"))
high_hm = Heatmap(as.matrix(high_corr), show_heatmap_legend = FALSE, 
	show_row_names = FALSE, show_column_names = FALSE, 
	show_column_dend = FALSE, col=col_fun,
	column_title = "Hot", column_title_gp = gpar(fontsize = 16))
low_hm = Heatmap(as.matrix(low_corr), name = "PCC", 
		show_row_names = FALSE, show_column_names = FALSE, col=col_fun,
		column_order = column_order(high_hm), row_order = row_order(high_hm),
		column_title = "Cold", column_title_gp = gpar(fontsize = 16))

hm_list = low_hm + high_hm
hmpng = paste(plot_dir, "merge_heatmap_ln_order_by_cluster_hot.png", sep = "") # Output: Fig2A
png(hmpng, res = 300, width = 12, heigh = 6, units = 'in')
draw(hm_list, main_heatmap = 2, row_dend_side = "right")
gar = dev.off()
cat("\t")
print(Sys.time() - st)


st = Sys.time()
cat("Ploting heatmaps (orderd by cold cohort)... \n")
low_hm = Heatmap(as.matrix(low_corr), name = "PCC", 
		show_row_names = FALSE, show_column_names = FALSE, 
		show_column_dend = FALSE, col = col_fun,
		column_title = "Cold", column_title_gp = gpar(fontsize = 16))

high_hm = Heatmap(as.matrix(high_corr), show_heatmap_legend = FALSE, 
	show_row_names = FALSE, show_column_names = FALSE, col = col_fun,
	column_order = column_order(low_hm), row_order = row_order(low_hm),
	column_title = "Hot", column_title_gp = gpar(fontsize = 16))
hm_list = low_hm + high_hm
hmpng = paste(plot_dir, "merge_heatmap_ln_order_by_cluster_cold.png", sep = "") # Output: Fig2A
png(hmpng, res = 300, width = 12, heigh = 6, units = 'in')
draw(hm_list, main_heatmap = 1, row_dend_side = "left")
gar = dev.off()
cat("\t")
print(Sys.time() - st)
