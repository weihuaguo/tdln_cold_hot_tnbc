rm(list = ls())

suppressMessages(library(ggplot2))
suppressMessages(library(ggridges))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(circlize))

download_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/ihc"# NOTE: Please change the ... into the directory you were saving the data

panel_dir <- paste(download_dir, "MC_IL4_mDC_Th2_Panel", sep = "/")
in_dir <- paste(panel_dir, "annotated_detections", sep = "/")
marker_loc <- c("IL4"="Cell", "CD117" = "Cytoplasm")
data_type <- "mean"
figres <- 300
print(names(marker_loc))

dtct_df <- read.table(paste(in_dir, "20240123_2250_30ROIs_detections.tsv", sep = "/"), check.names = F, sep = "\t", header = T)
dtct_df <- dtct_df[dtct_df$Name != "PathCellObject",]
info_df <- dtct_df[,c("Image", "Object ID", "Name", "Class", "Parent", "ROI", "Centroid X µm", "Centroid Y µm")]
data_df <- dtct_df[,str_detect(colnames(dtct_df), data_type)]
data_df <- data_df[,c("Cell: IL-4 (Opal 620) mean", "Cytoplasm: C-KIT (Opal 520) mean", "Cytoplasm: c-Kit (Opal 520) mean")]
data_df$`Cytoplasm: CD117 (Opal 520) mean` <- ifelse(is.na(data_df$`Cytoplasm: c-Kit (Opal 520) mean`), data_df$`Cytoplasm: C-KIT (Opal 520) mean`, data_df$`Cytoplasm: c-Kit (Opal 520) mean`)
clean_df <- cbind(info_df, data_df[,c("Cell: IL-4 (Opal 620) mean", "Cytoplasm: CD117 (Opal 520) mean")])
colnames(clean_df)[str_detect(colnames(clean_df), "mean")] <- c("IL4", "CD117")
clean_df$SID <- str_split_fixed(clean_df$Image, "_b", n = 2)[,1]
clean_df$ROIID <- str_c("ROI", str_split_fixed(clean_df$Image, "_ROI", n = 2)[,2])
clean_df$phenotype[clean_df$Class == "IL4"] <- "IL4+ non-MC"
clean_df$phenotype[clean_df$Class == "CD117"] <- "IL4- MC"
clean_df$phenotype[is.na(clean_df$phenotype)] <- "IL4+ MC"
write.csv(clean_df, paste(in_dir, "cleaned_used_cell_marker_mean_dataframe.csv", sep = "/"))

for (ifl in unique(clean_df$SID)) {
	cat(ifl, "\n")
	out_prf <- paste(ifl, data_type, "organized", sep = "_")
	df <- clean_df[clean_df$SID == ifl,]
	cts_df <- df %>%
		group_by(phenotype, Class) %>%
		summarize(n = n())
	bar_gg <- ggplot(cts_df, aes(x = phenotype, y = n, fill = Class)) +
		geom_bar(stat = "identity", color = "gray") +
		scale_y_continuous(trans='log10') +
		labs(y = "Cell number", x = "Annotated cell type", fill = "QuPath class\nbased on marker intensity") +
		theme_bw()
	ggsave(paste(in_dir, "/", out_prf, "_bar_cross_class.png", sep = ""), dpi = figres, width = 4.5, height = 6)

	merge_df <- df
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
		labs(x = "Average marker intensity (log2)", y = "Annotated cell type", title = ifl) +
		theme_bw() +
		theme(legend.position = "none")
	ggsave(paste(in_dir, "/", out_prf, "_ridge.png", sep = ""), dpi = figres, width = 6, height = 3)

	mean_df <- df[,names(marker_loc)]

	ha <- HeatmapAnnotation(Region = df$Parent, col = list(Region = c("BCellZone" = "purple", "TCellZone" = "green", "Sinus" = "goldenrod")))
	hm_obj <- Heatmap(t(scale(mean_df)),
		    name = "Z-score (mean)",
#		    col = color_scheme,
		    cluster_rows = TRUE,
		    cluster_columns = FALSE, 
		    column_names_gp = gpar(fontsize = 12),
#		    row_split = row_split_col,
		    row_title = "Marker",
		    column_split = df$phenotype,
		    column_order = order(df$Parent),
		    row_names_gp = gpar(fontsize = 12),
		    top_annotation = ha,
#		    row_order = row_orders,
		    show_column_names = F,
		    heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3.6, "cm")))
	png(paste(in_dir, "/", out_prf, "_precise_heatmap.png", sep = ""), res = figres, width = 9, height = 3, units = 'in')
	draw(hm_obj, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()
}
