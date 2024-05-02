# Compare and calculate correlation for each subset of genes
# Weihua Guo, Ph.D.
# 09/20/2021

rm(list = ls())
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(rstatix))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(ComplexHeatmap))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/nanostring"
smp_ann <- as.data.frame(read_excel(paste(data_dir, "/sample_annotation.xlsx", sep = "")))
raw_df <- as.data.frame(read_excel(paste(data_dir, "/merge_final_data.xlsx", sep = "")))
gene_ann <- as.data.frame(read_excel(paste(data_dir, "/immune_panel_cleaned_annotate.xlsx", sep = ""), 
				     sheet = "genome"))
de_df <- as.data.frame(read_excel(paste(data_dir, "/SD2_DE_NS_LN_list_240128.xlsx", sep = "")))
colnames(de_df)[1] <- "gene"

ppf <- paste(data_dir, "/clean_092021_ln_", sep = "")

rownames(smp_ann) <- smp_ann$pid
smp_ann$cohort <- ifelse(smp_ann$group == "high", "Hot", "Cold")
print(head(smp_ann))
rownames(raw_df) <- raw_df[,1]
raw_df[,1] <- NULL
print(raw_df[1:9,1:6])
rownames(gene_ann) <- gene_ann[,1]
gene_ann[,1] <- NULL
gene_ann <- gene_ann[gene_ann$`Gene Class` != "HK",]

clean_gene_ann <- gene_ann[,c("Annotation", "Common Name")]
print(head(clean_gene_ann))
clean_gene_ann$`Migration: chemokine` <- ifelse(str_detect(clean_gene_ann$Annotation, "Chemokines and receptors"), "+", "-")
clean_gene_ann$`Migration: cytokine` <- ifelse(str_detect(clean_gene_ann$Annotation, "Cytokines and receptors"), "+", "-")
clean_gene_ann$`Migration: leukocyte` <- ifelse(str_detect(clean_gene_ann$Annotation, "Leukocyte migration"), "+", "-")
gene_ann <- clean_gene_ann[,str_detect(colnames(clean_gene_ann), "Migration")]
gene_ann[is.na(gene_ann)] <- "-"
print(head(gene_ann))

ln_df <- raw_df[,raw_df["tissue",] == "LN"]
colnames(ln_df) <- ln_df["patient_id",]

avg_df <- as.data.frame(t(ln_df))
avg_df <- gather(avg_df, "gene", "expr", colnames(avg_df)[5:ncol(avg_df)])
avg_df$expr <- as.numeric(avg_df$expr)
avg_df <- avg_df %>%
	group_by(ds_subtype, gene) %>%
	summarize(mean = mean(expr, na.rm = T))
spr_avg_df <- spread(avg_df, "ds_subtype", "mean")
spr_avg_df <- merge(spr_avg_df, de_df, by = "gene", all.x = T)
print(head(spr_avg_df))

ln_df <- ln_df[5:nrow(ln_df),]
ln_df <- ln_df[,rownames(smp_ann)]
ln_df <- ln_df[rownames(gene_ann),]

print(head(ln_df))
write.csv(ln_df, paste(ppf, "tdln_log2_norm_expr.csv", sep = ""))
# Necessary for cor_test
ln_df <- read.csv(paste(ppf, "tdln_log2_norm_expr.csv", sep = ""), row.names=1)

use_df <- t(ln_df)
all_genes <- colnames(use_df)
use_df <- cbind(smp_ann, use_df)

hm_ann <- HeatmapAnnotation(Cohort = use_df$cohort, 
			    col = list(Cohort = c("Hot" = "firebrick", "Cold" = "dodgerblue")))

for (iss in colnames(gene_ann)) {
	cat(iss, "\n")
	tmp_genes <- rownames(gene_ann)[gene_ann[,iss] == "+"]
	print(length(tmp_genes))
	print(tmp_genes)
	print(head(gene_ann))
	tmp_avg_df <- merge(spr_avg_df, gene_ann, by.x = "gene", by.y = "row.names", all.x = T)
	tmp_avg_df$label <- ifelse(tmp_avg_df[,iss] == "+", tmp_avg_df$gene, NA)
	sct_gg <- ggplot(tmp_avg_df, aes(x = COLD, y = HOT, label = label)) +
		geom_point(aes(color = logFC)) +
		geom_label_repel() +
		scale_color_gradient2(midpoint = 0, low = "dodgerblue", mid = "white", high = "firebrick") +
		labs(title = iss, x = "Average expression (Cold TDLN)", y = "Average expression (Hot TDLN)") +
		theme_bw()
	ggsave(paste(ppf, iss, "_scatter_plot.png", sep = ""), sct_gg, dpi = 300, 
	       width = 4.5, height = 4.5, limitsize=FALSE)

	cat("\tDirectly compare...\n")
	tmp_df <- use_df[,tmp_genes]
	tmp_df <- cbind(smp_ann, tmp_df)
	tmp_gath_df <- gather(tmp_df, "gene", "expr", tmp_genes)
	tmp_gath_df$expr <- as.numeric(tmp_gath_df$expr)
	print(head(tmp_gath_df))
	comp_res <- tmp_gath_df %>%
		group_by(gene) %>%
		wilcox_test(expr ~ cohort, detailed = TRUE, p.adjust.method = "fdr") %>%
		adjust_pvalue(p.col = "p", method = "fdr") %>%
		add_significance() 
#	print(comp_res)
	write.csv(comp_res, paste(ppf, iss, "_wilcox_results.csv", sep = ""))

#	print(head(tmp_gath_df))
	row_num <- ceiling(length(tmp_genes)/6)
	print(row_num)
	box_gg <- ggplot(tmp_gath_df, aes(x = cohort, y = expr, color = cohort)) +
		geom_boxplot() +
		geom_point() +
		scale_color_manual(values=c("Cold" = "dodgerblue", "Hot" = "firebrick"))+
		scale_y_continuous(expand = c(0.1, 0.1, 0.2, 0.1)) +
		facet_wrap(~gene, ncol = 6, scale = "free") +
		stat_compare_means(label = 'p.format') +
		labs(x = "Cohort", y = "Expression\n(log2 transformed)", color = "Cohort", title = iss) +
		theme_bw() 
	ggsave(paste(ppf, iss, "_boxplot.png", sep = ""), box_gg, dpi = 300, 
	       width = 8.1, height = 2*row_num, limitsize=FALSE)

	cat("\tHeatmap...\n")
	tmp_df <- use_df[,tmp_genes]
#	for (i in 1:ncol(tmp_df)) {tmp_df[,i] <- as.numeric(tmp_df[,i])}
	hm_df <- t(scale(tmp_df))

	hm_obj <- Heatmap(hm_df,
			name = paste("Z-score\n(", iss, ")", sep = ""),
			cluster_rows = FALSE,
			cluster_columns = TRUE, 
			column_names_gp = gpar(fontsize = 9), 
			column_split = use_df$cohort,
			row_names_gp = gpar(fontsize = 9),
			top_annotation = hm_ann)
	png(paste(ppf, iss, "_expr_split_heatmap.png", sep = ""), 
	    res = 300, width = 6, height = length(tmp_genes)*0.25, units = 'in')
	print(draw(hm_obj, merge_legend = TRUE))
	gar <- dev.off()

	hm_obj <- Heatmap(hm_df,
			name = paste("Z-score\n(", iss, ")", sep = ""),
			cluster_rows = FALSE,
			cluster_columns = TRUE, 
			column_names_gp = gpar(fontsize = 9), 
			row_names_gp = gpar(fontsize = 9),
			top_annotation = hm_ann)
	png(paste(ppf, iss, "_expr_clst_heatmap.png", sep = ""), 
	    res = 300, width = 6, height = length(tmp_genes)*0.25, units = 'in')
	print(draw(hm_obj, merge_legend = TRUE))
	gar <- dev.off()

#	q(save = "no")
	cat("\tCalculate the correlation coefficient ")
	st <- Sys.time()
	cor_df <- use_df %>%
		group_by(cohort) %>%
		cor_test(vars = all_of(tmp_genes), 
			 vars2 = all_of(tmp_genes))
	cat("cost: ", Sys.time()-st, "\n")
	print(head(cor_df))
	write.csv(cor_df, paste(ppf, iss, "_correlation.csv", sep = ""))
	hot_cor_df <- cor_df[cor_df$cohort == "Hot",]
	hot_cor_df <- as.data.frame(spread(hot_cor_df[,c("var1", "var2", "cor")], "var2", "cor"))
	rownames(hot_cor_df) <- hot_cor_df[,1]
	hot_cor_df[,1] <- NULL
#	print(hot_cor_df)

	cold_cor_df <- cor_df[cor_df$cohort == "Cold",]
	cold_cor_df <- as.data.frame(spread(cold_cor_df[,c("var1", "var2", "cor")], "var2", "cor"))
	rownames(cold_cor_df) <- cold_cor_df[,1]
	cold_cor_df[,1] <- NULL
#	print(cold_cor_df)

	hot_hm <- Heatmap(hot_cor_df, name = "PCC", show_column_dend=FALSE, 
			  column_title = "Hot", column_title_gp = gpar(fontsize = 16))
	cold_ho_hm <- Heatmap(cold_cor_df, name = "PCC", 
			      column_order = column_order(hot_hm), row_order = row_order(hot_hm),
			      column_title = "Cold (ordered by hot)",column_title_gp = gpar(fontsize = 16))
	ho_list <- cold_ho_hm+hot_hm
	hmpng <- paste(ppf, iss, "_hot_ordered_corr_hm.png", sep = "")
	png(hmpng, res = 300, width = length(tmp_genes)*0.5, heigh = length(tmp_genes)*0.25, units = 'in')
	print(draw(ho_list, main_heatmap = 2, row_dend_side = "right"))
	gar <- dev.off()

	cold_hm <- Heatmap(cold_cor_df, name = "PCC", show_column_dend=FALSE, 
			  column_title = "Cold", column_title_gp = gpar(fontsize = 16))
	hot_co_hm <- Heatmap(hot_cor_df, name = "PCC", 
			      column_order = column_order(cold_hm), row_order = row_order(cold_hm),
			      column_title = "Hot (ordered by cold)",column_title_gp = gpar(fontsize = 16))
	co_list <- cold_hm+hot_co_hm
	hmpng <- paste(ppf, iss, "_cold_ordered_corr_hm.png", sep = "")
	png(hmpng, res = 300, width = length(tmp_genes)*0.5, heigh = length(tmp_genes)*0.25, units = 'in')
	print(draw(co_list, main_heatmap = 1, row_dend_side = "right"))
	gar <- dev.off()
}
