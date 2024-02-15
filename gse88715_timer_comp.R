rm(list = ls())

suppressMessages(library(dplyr))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(rstatix))
suppressMessages(library(stringr))
suppressMessages(library(ggrepel))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tnbc_ln_cold_hot_ns/nanostring_results/park_jci"
data_file <- "estimation_matrix.csv"
group_file <- "gsm_pid_immg_jci.xlsx"

df <- read.csv(paste(data_dir, data_file, sep = "/"), header = T, row.names = 1)
grp_df <- as.data.frame(read_excel(paste(data_dir, group_file, sep = "/")))
df$ct_tl <- rownames(df)
df$ct <- str_split_fixed(rownames(df), "_", n = 2)[,1]
df$tl <- str_split_fixed(rownames(df), "_", n = 2)[,2]

gath_df <- gather(df, "gsm_id", "value", colnames(df)[str_detect(colnames(df), "GSM")])
gath_df <- merge(gath_df, grp_df, by = "gsm_id", all.x = T)
print(head(gath_df))

for (itl in unique(gath_df$tl)) {
	cat(itl, "\n")
	tmp_gath_df <- gath_df[gath_df$tl == itl,]
	if (nrow(tmp_gath_df) > 0) {
	tmp_gath_df$immgroup <- factor(tmp_gath_df$immgroup, levels = c("ID", "MR", "SR", "FI"))
#	print(tmp_sgnt)
	print(head(tmp_gath_df))

	cat("\tTwo groups\n")
	tmp_gath_df$group <- ifelse(tmp_gath_df$immgroup == "ID" | tmp_gath_df$immgroup == "MR", "Cold (ID+MR)", "Hot (SR+FI)")
	tmp_gath_df$cohort <- ifelse(tmp_gath_df$immgroup == "ID" | tmp_gath_df$immgroup == "MR", "Cold", "Hot")

	stat_df <- tmp_gath_df %>% 
		group_by(ct, tissue, group) %>%
		summarise(n = n(),
			  mean = mean(value, na.rm = T),
			  median = median(value, na.rm = T),
			  sd = sd(value, na.rm = T),
			  min = min(value, na.rm = T),
			  max = max(value, na.rm = T),
			  se = sd/sqrt(n),
			  mean_se_upper = mean+se,
			  mean_se_lower = mean-se)
	write.csv(stat_df, paste(data_dir, "/timer20_2g_", itl, "_descriptive_statistics.csv", sep = ""))
	sd_zero <- str_c(stat_df$ct, stat_df$tissue)[stat_df$sd != 0]

	sum_timer_df <- tmp_gath_df %>%
		group_by(ct, cohort, tissue) %>%
		summarize(mean = mean(value, na.rm = T)) %>%
		spread(cohort, mean)
	sum_timer_df$avg_log2FC <- log2(sum_timer_df$Cold/sum_timer_df$Hot)
#	print(head(sum_timer_df))
	boxgg <- ggplot(tmp_gath_df, aes(x = ct, y = value, color = group)) +
		geom_boxplot() +
		geom_point(position = position_dodge(width = 0.75), alpha = 0.6, size = 1) +
		stat_compare_means(aes(group = group), label = "p.format") +
		facet_wrap(~tissue, nrow = 1) +
		scale_color_manual(values = c("Cold (ID+MR)" = "dodgerblue", "Hot (SR+FI)" = "firebrick")) +
		labs(y = itl) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggsave(paste(data_dir, "/timer20_2g_", itl, "_box.png", sep = ""), boxgg, dpi = 300, 
	       width = 2+1.8*length(unique(tmp_gath_df$ct)), height = 6, limitsize=F)

	tmp_gath_df$sd_zero_id <- str_c(tmp_gath_df$ct, tmp_gath_df$tissue)
	pw_res <- tmp_gath_df %>% 
		filter(sd_zero_id %in% sd_zero) %>%
		group_by(ct, tissue) %>%
		pairwise_wilcox_test(value ~ group, detailed = T) %>%
		add_significance('p') %>%
		adjust_pvalue('p') %>%
		add_significance('p.adj')

	pw_res <- merge(pw_res, sum_timer_df, by = c("ct", "tissue"), all.x = T)
	write.csv(pw_res, paste(data_dir, "/timer20_2g_", itl, "_gene_expr_pairwise_wilcox.csv", sep = ""))

	pw_res$logp <- -log10(pw_res$p)
	pw_res$logpadj <- -log10(pw_res$p.adj)
	vlcgg <- ggplot(pw_res, aes(x = avg_log2FC, y = logp)) +
		geom_point() +
		geom_text_repel(aes(label = ct)) +
		geom_vline(xintercept = c(-1,1), linetype = "dashed") +
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), linetype = "dashed") +
		facet_wrap(~tissue, nrow = 1, scales = "free") +
		labs(x = "Hot\t\t\tlog2 fold change (average)\t\t\tCold", y = "-log10 p-value", title = itl) +
		theme_classic()
	ggsave(paste(data_dir, "/timer20_2g_", itl, "_timer_comparison_volcano_p.png", sep = ""), dpi = 300, width = 10.8, height = 6, limitsize = F)
	vlcgg <- ggplot(pw_res, aes(x = avg_log2FC, y = logpadj)) +
		geom_point() +
		geom_text_repel(aes(label = ct)) +
		geom_vline(xintercept = c(-1,1), linetype = "dashed") +
		geom_vline(xintercept = 0) +
		geom_hline(yintercept = c(-log10(0.05), -log10(0.01)), linetype = "dashed") +
		facet_wrap(~tissue, nrow = 1, scales = "free") +
		labs(x = "Hot\t\t\tlog2 fold change (average)\t\t\tCold", y = "-log10 adjusted p-value", title = itl) +
		theme_classic()
	ggsave(paste(data_dir, "/timer20_2g_", itl, "_timer_comparison_volcano_padj.png", sep = ""), dpi = 300, width = 10.8, height = 6, limitsize = F)

	tmp_df <- df[df$tl == itl,]
	rownames(tmp_df) <- tmp_df$ct
	group_df <- grp_df[grp_df$tissue == "Stroma",]
	group_df$cohort <- ifelse(group_df$immgroup == "ID" | group_df$immgroup == "MR", "Cold", "Hot")
	hm_top <- pw_res[pw_res$tissue == "Stroma",]

	hm_expr_df <- tmp_df[hm_top$ct, group_df$gsm_id]
#	print(head(hm_expr_df))

	padjColFun <- c("ns" = "grey", "*"="darkgoldenrod1", "**"="darkgoldenrod2", "***"="darkgoldenrod3", "****"="darkgoldenrod4" )
	gene_ann <- rowAnnotation(log2FC = anno_barplot(hm_top$avg_log2FC, baseline=0, ylim = c(-max(abs(hm_top$avg_log2FC)), max(abs(hm_top$avg_log2FC)))),
				  AdjP = hm_top$p.adj.signif,
				      col = list(AdjP = padjColFun),
				      annotation_legend_param = list(AdjP = list(title = "Adjusted P-value", legend_height = unit(2, "cm"), nrow = 1,
										 legend_direction = "horizontal", gp=gpar(fontsize = 14)))
	)
	smp_ann <- HeatmapAnnotation(Cohort = group_df$cohort, col = list(Cohort = c("Hot" = "firebrick", "Cold" = "dodgerblue")),
				      annotation_legend_param = list(Cohort = list(title = "Cohort", legend_height = unit(2, "cm"), nrow = 1,
										 legend_direction = "horizontal", gp=gpar(fontsize = 14)))

	)

	hm <- Heatmap(t(scale(t(hm_expr_df))),
			 name = "Z-score",
			 cluster_rows = TRUE,
			 cluster_columns = TRUE,
			 column_names_gp = gpar(fontsize = 9), 
			 row_names_gp = gpar(fontsize = 12),
			 left_annotation = gene_ann,
			 top_annotation = smp_ann,
			 show_column_names = F,
			 column_title = paste("Stroma", itl),
			 heatmap_legend_param = list(direction = "horizontal", legend_height = unit(3.6, "cm")))
	png(paste(data_dir, "/timer20_2g_", itl, "_stroma_heatmap.png", sep = ""), res=300, width=6, height=nrow(hm_expr_df)/2.7, units='in')
	draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
	gar <- dev.off()

	tmp_df <- df[df$tl == itl,]
	rownames(tmp_df) <- tmp_df$ct
	group_df <- grp_df[grp_df$tissue == "Epithelium",]
	group_df$cohort <- ifelse(group_df$immgroup == "ID" | group_df$immgroup == "MR", "Cold", "Hot")
	hm_top <- pw_res[pw_res$tissue == "Epithelium",]

	hm_expr_df <- tmp_df[hm_top$ct, group_df$gsm_id]
#	print(head(hm_expr_df))
#	print(nrow(hm_expr_df))
	hm_top$avg_log2FC[hm_top$avg_log2FC == Inf] <- max(hm_top$avg_log2FC[hm_top$avg_log2FC != Inf])*1.1
	print(head(hm_top))

	padjColFun <- c("ns" = "grey", "*"="darkgoldenrod1", "**"="darkgoldenrod2", "***"="darkgoldenrod3", "****"="darkgoldenrod4" )
	gene_ann <- rowAnnotation(log2FC = anno_barplot(hm_top$avg_log2FC, baseline=0, ylim = c(-max(abs(hm_top$avg_log2FC)), max(abs(hm_top$avg_log2FC)))),
				  AdjP = hm_top$p.adj.signif,
				      col = list(AdjP = padjColFun),
				      annotation_legend_param = list(AdjP = list(title = "Adjusted P-value", legend_height = unit(2, "cm"), 
										 legend_direction = "vertical", gp=gpar(fontsize = 14))
								     ),
				      annotation_name_gp= gpar(fontsize = 10)
	)
	smp_ann <- HeatmapAnnotation(Cohort = group_df$cohort, col = list(Cohort = c("Hot" = "firebrick", "Cold" = "dodgerblue")), annotation_legend_param = list(gp=gpar(fontsize = 14)),
				     annotation_name_gp= gpar(fontsize = 10)
	)

	hm <- Heatmap(t(scale(t(hm_expr_df))),
			 name = "Z-score",
			 cluster_rows = TRUE,
			 cluster_columns = TRUE,
			 column_names_gp = gpar(fontsize = 9), 
			 row_names_gp = gpar(fontsize = 7.5),
			 left_annotation = gene_ann,
			 top_annotation = smp_ann,
			 show_column_names = F,
			 column_title = paste("Epithelium", itl),
			 heatmap_legend_param = list(direction = "vertical", legend_height = unit(3.6, "cm")))
	png(paste(data_dir, "/timer20_2g_", itl, "_epithelium_heatmap.png", sep = ""), res=300, width=9, height=nrow(hm_expr_df)/2.7, units='in')
	draw(hm, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
	gar <- dev.off()
	}
}
