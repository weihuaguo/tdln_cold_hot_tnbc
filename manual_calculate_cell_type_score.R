rm(list = ls())
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(rstatix))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(ComplexHeatmap))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/nanostring"
input_file <- "ImmunePanel_ALL_normalized_data_LN_FINAL.xlsx"
smp_ann <- as.data.frame(read_excel(paste(data_dir, "/sample_annotation.xlsx", sep = "")))
raw_df <- as.data.frame(read_excel(paste(data_dir, input_file, sep = "/")))
gene_ann <- as.data.frame(read_excel(paste(data_dir, "/immune_panel_cleaned_annotate_v2.xlsx", sep = ""), 
				     sheet = "genome"))
#gene_ann <- as.data.frame(read_excel(paste(data_dir, "/profile_panel_cleaned_annotate_v1.xlsx", sep = ""), 
#				     sheet = "genome"))

exp_id <- "clean_230809_ln_immune_profile"
print(raw_df[1:9,1:6])
pf <- paste(data_dir, "/", exp_id, "_", sep = "")

rownames(raw_df) <- raw_df$pid
raw_df$pid <- NULL

gene_ann <- gene_ann[!is.na(gene_ann$`Cell Type`),]
cts_df <- as.data.frame(matrix(ncol = length(unique(gene_ann$`Cell Type`)), nrow = nrow(raw_df)))
rownames(cts_df) <- rownames(raw_df)
colnames(cts_df) <- unique(gene_ann$`Cell Type`)
for (ic in unique(gene_ann$`Cell Type`)) {
	cat(ic, "\n")
	tmp_genes <- gene_ann$`HUGO Name`[gene_ann$`Cell Type` == ic]
	print(tmp_genes)
	ava_genes <- intersect(tmp_genes, colnames(raw_df))
	print(ava_genes)
	diff_genes <- setdiff(tmp_genes, colnames(raw_df))
	print(diff_genes)
	if (length(ava_genes) > 1) {
	cts_df[,ic] <- rowSums(raw_df[,ava_genes])/length(ava_genes)
	} else if (length(ava_genes) == 1) {
		cts_df[,ic] <- raw_df[,ava_genes]
	}
}
merge_cts_df <- merge(cts_df, smp_ann, by.x = "row.names", by.y = "pid", all.x = T)
write.csv(merge_cts_df, paste(pf, "merge_manual_cell_type_score.csv", sep = ""))
