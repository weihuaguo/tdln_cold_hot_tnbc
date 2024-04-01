rm(list = ls())

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/public_data/gse139568"
smp_file <- "sample_gse139568.xlsx"
exp_file <- "expr_gse139568.xlsx"
gn_file <- "genome_gse139568.xlsx"
ann_file <- "immune_panel_cleaned_annotate_v2.xlsx"
save_pf <- paste(data_dir, "/median_cell_type_", sep = "")

expr_df <- as.data.frame(read_excel(paste(data_dir, exp_file, sep = "/")))
colnames(expr_df)[1] <- "ID"
print(dim(expr_df))
print(expr_df[1:9,1:6])

gn_df <- as.data.frame(read_excel(paste(data_dir, gn_file, sep = "/")))
print(head(gn_df))

smp_df <- as.data.frame(read_excel(paste(data_dir, smp_file, sep = "/")))
print(head(smp_df))

ann_df <- as.data.frame(read_excel(paste(data_dir, ann_file, sep = "/")))
ann_df <- ann_df[!is.na(ann_df$`Cell Type`),]
ann_df$Name <- str_c(str_sub(ann_df$`HUGO Name`,1,1), tolower(str_sub(ann_df$`HUGO Name`, 2,-1)))
print(head(ann_df))

c <- 0

for (ic in unique(ann_df$`Cell Type`)) {
	cat(ic, "\n")
	tmp_ann <- ann_df[ann_df$`Cell Type` == ic,]
	print(tmp_ann)
	tmp_gn <- gn_df[gn_df$GENE_SYMBOL %in% tmp_ann$Name,]
	tmp_ann <- merge(tmp_ann, gn_df, by.x = "Name", by.y = "GENE_SYMBOL", all.x = T)

	print(dim(tmp_gn))
	print(dim(tmp_ann))

	tmp_exp <- merge(tmp_ann, expr_df, by = "ID", all.x = T)
	write.csv(tmp_exp, paste(save_pf, ic, "_everything.csv", sep = ""))
	print(dim(tmp_exp))

	tmp_gath <- gather(tmp_exp, "SID", "expr", colnames(tmp_exp)[str_detect(colnames(tmp_exp), "GSM")])
	tmp_gath <- merge(tmp_gath, smp_df, by = "SID", all.x = T)
	tmp_gath$GID <- str_c(tmp_gath$Name, ">", tmp_gath$ID)
	print(head(tmp_gath))
	vln_gg <- ggplot(tmp_gath, aes(x = GID, y = expr, color = group)) +
		geom_boxplot() +
		geom_point(position = position_dodge(width = 0.7)) +
		stat_compare_means(aes(group = group), label = "p.format") +
		labs(x = "Gene>ID", y = "Expression", title = ic, color = "Group") +
		theme_classic() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	ggsave(paste(save_pf, ic, "_boxplot.png", sep = ""), dpi = 300, width = 2+0.7*nrow(tmp_ann), height = 6, limitsize=F)
	tmp_score <- tmp_gath %>%
		group_by(group, SID, `Cell Type`) %>%
		summarize(score = mean(expr, na.rm = T))
	print(head(tmp_score))
	if (c == 0) {
		score_df <- tmp_score
	} else {
		score_df <- rbind(score_df, tmp_score)
	}
	c <- c+1
}
write.csv(score_df, paste(save_pf, "ALL_scores.csv", sep = ""))

vln_gg <- ggplot(score_df, aes(x = `Cell Type`, y = score, color = group)) +
	geom_boxplot() +
	geom_point(position = position_dodge(width = 0.7)) +
	stat_compare_means(aes(group = group), label = "p.format") +
	labs(y = "Score", title = ic, color = "Group") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(save_pf,"ALL_boxplot.png", sep = ""), dpi = 300, width = 2+0.7*length(unique(ann_df$`Cell Type`)), height = 6, limitsize=F)

