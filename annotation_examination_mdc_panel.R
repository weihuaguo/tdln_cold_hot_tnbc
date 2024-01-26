rm(list = ls())

suppressMessages(library(ggplot2))
suppressMessages(library(ggridges))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(circlize))

download_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/ihc"# NOTE: Please change the ... into the directory you were saving the data

panel_dir <- paste(download_dir, "mDC_Th_Panel", sep = "/")
in_dir <- paste(panel_dir, "annotated_detections", sep = "/")
marker_loc <- c("CD3"="Cytoplasm", "CD4" = "Cytoplasm", "GATA3" = "Nucleus", "DC-LAMP" = "Cell", "Tbet" = "Nucleus")
data_type <- "mean"
figres <- 300
print(names(marker_loc))

all_files <- list.files(in_dir, pattern = "extracted_detection.csv")
print(all_files)
for (ifl in all_files) {
	cat(ifl, "\n")
	sid <- str_split_fixed(ifl, "_ex", n = 2)[1]
	out_prf <- paste(sid, data_type, "organized", sep = "_")
	df <- read.csv(paste(in_dir, ifl, sep = "/"), check.names = F, row.names = 1)
	df$phenotype <- ifelse(df$phenotype == "CD8+", "CD3+CD4-", df$phenotype)
	print(head(df))
	print(colnames(df))

	mean_df <- as.data.frame(matrix(ncol = length(names(marker_loc)), nrow = nrow(df)))
	colnames(mean_df) <- names(marker_loc)
	rownames(mean_df) <- rownames(df)
	for (im in names(marker_loc)) {
		cat("\t", im, "\n")
		tmp_col <- paste(marker_loc[im], im, sep = ": ")
		df_col <- colnames(df)[str_detect(colnames(df), "mean") & str_detect(colnames(df), tmp_col)]
		mean_df[,im] <- df[,df_col]
	}
	write.csv(mean_df, paste(in_dir, "/", out_prf, "_extracted_mean_vis.csv", sep = ""))
	print(head(df))

	cts_df <- df %>%
		group_by(phenotype, Class) %>%
		summarize(n = n())
	bar_gg <- ggplot(cts_df, aes(x = phenotype, y = n, fill = Class)) +
		geom_bar(stat = "identity", color = "gray") +
		scale_y_continuous(trans='log10') +
		labs(y = "Cell number", x = "Annotated cell type", fill = "QuPath class\nbased on marker intensity") +
		theme_bw()
	ggsave(paste(in_dir, "/", out_prf, "_bar_cross_class.png", sep = ""), dpi = figres, width = 4.5, height = 6)

	merge_df <- cbind(mean_df, df[,c("phenotype", "Class")])
	gath_df <- gather(merge_df, "marker", "mean", names(marker_loc))
	gath_df$log_mean <- log2(gath_df$mean + 1)
	print(head(gath_df))

	rdg_gg <- ggplot(gath_df, aes(x = log_mean, y = phenotype)) +
		geom_density_ridges(jittered_points = TRUE, quantile_lines = TRUE, scale = 0.9, alpha = 0.7, 
				    vline_size = 1, vline_color = "red",
				    position = position_points_jitter(width = 0.05, height = 0, adjust_vlines = TRUE),
				    point_shape = '|', point_size = 0.9, point_alpha = 0.05,
				    aes(fill = phenotype)
		) +
		facet_wrap(~marker, nrow = 1, scales = "free_x") +
		labs(x = "Average marker intensity (log2)", y = "Annotated cell type", title = sid) +
		theme_bw() +
		theme(legend.position = "none")
	ggsave(paste(in_dir, "/", out_prf, "_ridge.png", sep = ""), dpi = figres, width = 9, height = 4)

	hm_obj <- Heatmap(t(scale(mean_df)),
		    name = "Z-score (mean)",
#		    col = color_scheme,
		    cluster_rows = TRUE,
		    cluster_columns = FALSE, 
		    column_names_gp = gpar(fontsize = 12),
#		    row_split = row_split_col,
		    row_title = paste("Marker\n(", sid, ")", sep = ""),
		    column_split = df$phenotype,
		    column_order = order(df$Class),
		    row_names_gp = gpar(fontsize = 12),
#		    top_annotation = col_hm_ann,
#		    row_order = row_orders,
		    show_column_names = F,
		    heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(in_dir, "/", out_prf, "_precise_heatmap.png", sep = ""), res = figres, width = 12, height = 3, units = 'in')
	draw(hm_obj, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()
}
