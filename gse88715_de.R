rm(list = ls())

res_out = function(fit, gann, coef = "all", save_dir = "", fitname = "", plot_dir = "") {
	cat("\tStart to output results from limma...\n")
	cont_name = colnames(fit$contrasts)
	if (coef != "all") {
		cat("\t  Extracting results for contrasts ", cont_name[coef], "\n")
		restable = topTable(fit, coef = coef, number=nrow(fit), adjust.method="fdr")
		if (save_dir == "") {
			cat("\t   WARNING: No files are saved!!!\n")
		} else {
			cat("\t   Result table will be saved to: \n")
			csv_file = paste(save_dir, fitname, cont_name[coef], ".csv", sep = "_")
			cat("\t     ", csv_file, "\n")
			print_res = merge(gene_annot, restable, by = 0)
			write.csv(print_res, file = csv_file)
		}
	} else {
		for (i in 1:ncol(fit$contrasts)) {
			cat("\t  Extracting results for contrasts ", cont_name[i], "\n")
			restable = topTable(fit, coef = i, number=nrow(fit), adjust.method="fdr")
			if (save_dir == "") {
				cat("\t   WARNING: No files are saved!!!\n")
			} else {
				cat("\t   Result table will be saved to: \n")
				csv_file = paste(save_dir, fitname, cont_name[i], ".csv", sep = "_")
				cat("\t     ", csv_file, "\n")
				print_res = merge(gene_annot, restable, by = 0)
				write.csv(print_res, file = csv_file)
				if (plot_dir == "") {
					cat("\t     NO Volcano plots!!!\n")
				} else {
					cat("\t   Enhanced Volcano plots will be save to: \n")
					tiff_file = paste(save_dir, fitname, cont_name[i], ".tiff", sep = "_")
					cat("\t     ", tiff_file, "\n")

					evplot = EnhancedVolcano(print_res, 
						x = 'logFC', y = 'adj.P.Val',
						lab = print_res$GENE_SYMBOL,
						pCutoff = 0.10, FCcutoff = 1.5, 
						xlim = c(-8, 8), # ylim = c(0, 6),
						title = paste(fitname, cont_name[i], sep = " "),
						transcriptPointSize = 1.5,
						transcriptLabSize = 3.,
						xlab = bquote(~Log[2]~ 'fold change'),
						ylab = bquote(~-Log[10]~adjusted~italic(P)),
						cutoffLineWidth = 0.8,
						cutoffLineCol = 'black',
						legend=c('NS','Log (base 2) fold-change','Adjusted P value\n (FDR)',
							'Significantly \ndifferentially \nexpressed genes'),
						legendPosition = 'right')
					tiff(tiff_file, res = 120, width = 9, heigh = 6, units = 'in')
					print(evplot)
					gar = dev.off()
				}
			}
		}
	}
}

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(stringr))
suppressMessages(library(limma))
suppressMessages(library(circlize))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tnbc_ln_cold_hot_ns/nanostring_results/park_jci"
data_file <- "jci_expr_with_annotation_v4.xlsx"
group_file <- "gsm_pid_immg_jci.xlsx"

raw_expr <- as.data.frame(read_excel(paste(data_dir, data_file, sep = "/")))
gene_annot <- raw_expr[,1:3]
rownames(gene_annot) <- gene_annot[,1]

clean_expr <- raw_expr[,-c(2,3)]
rownames(clean_expr) <- clean_expr[,1]
clean_expr <- clean_expr[,-c(1)]
cat("Cleaned expression matrix: ", dim(clean_expr), "\n")

### Read the grouping of each sample
sgroup <- read_excel(paste(data_dir, group_file, sep = "/"))
sgroup <- as.data.frame(sgroup)
rownames(sgroup) <- sgroup[,1]
sgroup <- sgroup[,-1]

if (all(rownames(sgroup) == colnames(clean_expr))) {
	cat("Annotation matrix and expression matrix is matched!\n")
} else {stop("The annotation and expression matrix are NOT matched!!!")}

sgroup$group <- ifelse(sgroup$immgroup == "ID" | sgroup$immgroup == "MR", "Cold", "Hot")
sgroup$comb_fact <- paste(sgroup$tissue, sgroup$group, sep =".")
tann <- factor(sgroup$tissue)
iann <- factor(sgroup$group)
comb_ann <- factor(sgroup$comb_fact)

design <- model.matrix(~0+comb_ann)
colnames(design) <- levels(comb_ann)
fit <- lmFit(clean_expr, design)

# For Cancer Islands (CI)
cat("Compare different immune groups in cancer island...\n")
cicont.mat <- makeContrasts(
	cvsh = Epithelium.Cold-Epithelium.Hot,
	levels = design)
cifit2 <- contrasts.fit(fit, cicont.mat)
cifit2 <- eBayes(cifit2)

restable <- topTable(cifit2, number=nrow(cifit2), adjust.method="fdr")
restable <- merge(restable, gene_annot, by = "row.names", all.x = T)
print(head(restable))
write.csv(restable, paste(data_dir, "/limma_2g_ci_result.csv", sep = ""))

evplot = EnhancedVolcano(restable, 
	x = 'logFC', y = 'adj.P.Val',
	lab = restable$GENE_SYMBOL,
	pCutoff = 0.10, FCcutoff = 1.0, 
	title = "Within Cancer Islands (Cold [ID+MR] vs Hot [SR+FI])",
	xlab = bquote(~Log[2]~ 'fold change'),
	ylab = bquote(~-Log[10]~adjusted~italic(P)),
	legendPosition = 'right')
png(paste(data_dir, "/limma_2g_ci_ev_padj.png", sep = ""), res = 300, width = 9, heigh = 6, units = 'in')
print(evplot)
gar <- dev.off()

evplot = EnhancedVolcano(restable, 
	x = 'logFC', y = 'P.Value',
	lab = restable$GENE_SYMBOL,
	pCutoff = 0.10, FCcutoff = 1.0, 
	title = "Within Cancer Islands (Cold [ID+MR] vs Hot [SR+FI])",
	xlab = bquote(~Log[2]~ 'fold change'),
	ylab = bquote(~-Log[10]~italic(P)),
	legendPosition = 'right')
png(paste(data_dir, "/limma_2g_ci_ev_p.png", sep = ""), res = 300, width = 9, heigh = 6, units = 'in')
print(evplot)
gar <- dev.off()

cat("Heatmaps\n")
up_top <- restable %>% filter(adj.P.Val <= 0.10, logFC>=1) %>% top_n(wt = logFC, n = 20)
print(up_top)
down_top <- restable %>% filter(adj.P.Val <= 0.10, logFC<=-1) %>% top_n(wt = logFC, n = -20)
print(down_top)

hm_top <- rbind(up_top, down_top)
hm_expr_df <- clean_expr[hm_top$Row.names,]

print(hm_expr_df[1:9,1:6])
print(clean_expr[1:9,1:6])

fcColFun <- colorRamp2(breaks = c(min(hm_top$logFC), 0, max(hm_top$logFC)), colors = c("orange", "white", "forestgreen"))
gene_ann <- rowAnnotation(log2FC = hm_top$logFC,
			      col = list(log2FC = fcColFun),
			      annotation_legend_param = list(log2FC = list(title = "log2FC", legend_height = unit(2, "cm"), legend_direction = "vertical", gp=gpar(fontsize = 14))),
			      annotation_name_gp= gpar(fontsize = 10)
)

hm_group <- sgroup[sgroup$tissue == "Epithelium",]
print(head(sgroup))
smp_ann <- HeatmapAnnotation(Cohort = hm_group$group, 
			     col = list(Cohort = c("Hot" = "firebrick", "Cold" = "dodgerblue")),
			     annotation_legend_param = list(gp=gpar(fontsize = 14)),
			     annotation_name_gp= gpar(fontsize = 10)
)

rownames(hm_expr_df) <- hm_top$GENE_SYMBOL
hm_expr_df <- hm_expr_df[,rownames(hm_group)]
hm <- Heatmap(t(scale(t(hm_expr_df))),
                 name = "Z-score",
                 cluster_rows = TRUE,
                 cluster_columns = TRUE,
		 show_column_names = FALSE,
                 column_names_gp = gpar(fontsize = 9), 
                 row_names_gp = gpar(fontsize = 7.5),
                 left_annotation = gene_ann,
                 top_annotation = smp_ann,
                 heatmap_legend_param = list(direction = "vertical", legend_height = unit(3.6, "cm")))
png(paste(data_dir, "/limma_2g_ci_top20_heatmap.png", sep = ""), res=300, width=4.5, height=4.5, units='in')
draw(hm, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
gar <- dev.off()


cat("Compare different immune groups in stroma...\n")
stcont.mat <- makeContrasts(
	cvsh = Stroma.Cold-Stroma.Hot,
	levels = design)
stfit2 <- contrasts.fit(fit, stcont.mat)
stfit2 <- eBayes(stfit2)

restable <- topTable(stfit2, number=nrow(stfit2), adjust.method="fdr")
restable <- merge(restable, gene_annot, by = "row.names", all.x = T)
print(head(restable))
write.csv(restable, paste(data_dir, "/limma_2g_stroma_result.csv", sep = ""))

restable$logFC <- -restable$logFC
evplot = EnhancedVolcano(restable, 
	x = 'logFC', y = 'adj.P.Val',
	lab = restable$GENE_SYMBOL,
	pCutoff = 0.10, FCcutoff = 1.0, 
	xlim = c(-4, 3), ylim = c(0, 4.5),
	title = "Within Stroma (Cold [ID+MR] vs Hot [SR+FI])",
	xlab = bquote(~Log[2]~ 'fold change'),
	ylab = bquote(~-Log[10]~adjusted~italic(P)),
	drawConnectors = TRUE,
	legendPosition = 'right')
png(paste(data_dir, "/limma_2g_stroma_ev_padj.png", sep = ""), res = 300, width = 9, heigh = 6, units = 'in')
print(evplot)
gar <- dev.off()

cat("Heatmaps\n")
up_top <- restable %>% filter(adj.P.Val <= 0.10, logFC>=1) %>% top_n(wt = logFC, n = 20)
print(up_top)
down_top <- restable %>% filter(adj.P.Val <= 0.10, logFC<=-1) %>% top_n(wt = logFC, n = -20)
print(down_top)

hm_top <- rbind(up_top, down_top)
hm_expr_df <- clean_expr[hm_top$Row.names,]

print(hm_expr_df[1:9,1:6])
print(clean_expr[1:9,1:6])

fcColFun <- colorRamp2(breaks = c(min(hm_top$logFC), 0, max(hm_top$logFC)), colors = c("orange", "white", "forestgreen"))
gene_ann <- rowAnnotation(log2FC = hm_top$logFC,
			      col = list(log2FC = fcColFun),
			      annotation_legend_param = list(log2FC = list(title = "log2FC", legend_height = unit(2, "cm"), legend_direction = "vertical", gp=gpar(fontsize = 14))),
			      annotation_name_gp= gpar(fontsize = 10)
)

hm_group <- sgroup[sgroup$tissue == "Stroma",]
print(head(sgroup))
smp_ann <- HeatmapAnnotation(Cohort = hm_group$group, 
			     col = list(Cohort = c("Hot" = "firebrick", "Cold" = "dodgerblue")),
			     annotation_legend_param = list(gp=gpar(fontsize = 14)),
			     annotation_name_gp= gpar(fontsize = 10)
)

rownames(hm_expr_df) <- hm_top$GENE_SYMBOL
hm_expr_df <- hm_expr_df[,rownames(hm_group)]
hm <- Heatmap(t(scale(t(hm_expr_df))),
                 name = "Z-score",
                 cluster_rows = TRUE,
                 cluster_columns = TRUE,
		 show_column_names = FALSE,
                 column_names_gp = gpar(fontsize = 9), 
                 row_names_gp = gpar(fontsize = 7.5),
                 left_annotation = gene_ann,
                 top_annotation = smp_ann,
                 heatmap_legend_param = list(direction = "vertical", legend_height = unit(3.6, "cm")))
png(paste(data_dir, "/limma_2g_stroma_top20_heatmap.png", sep = ""), res=300, width=4.5, height=4.5, units='in')
draw(hm, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")
gar <- dev.off()

evplot = EnhancedVolcano(restable, 
	x = 'logFC', y = 'P.Value',
	lab = restable$GENE_SYMBOL,
	pCutoff = 0.10, FCcutoff = 1.0, 
	title = "Within Stroma (Cold [ID+MR] vs Hot [SR+FI])",
	xlab = bquote(~Log[2]~ 'fold change'),
	ylab = bquote(~-Log[10]~italic(P)),
	legendPosition = 'right')
png(paste(data_dir, "/limma_2g_stroma_ev_p.png", sep = ""), res = 300, width = 9, heigh = 6, units = 'in')
print(evplot)
gar <- dev.off()
