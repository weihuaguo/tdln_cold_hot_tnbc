
rm(list = ls())

suppressMessages(library(limma))
suppressMessages(library(dplyr))
suppressMessages(library(readxl))
suppressMessages(library(xlsx))
suppressMessages(library(EnhancedVolcano))

download_dir <- "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir <- paste(download_dir, "input_file", sep = "/")
out_dir <- paste(download_dir, "FigS1", sep = "/")
expr_xlsx <- paste(in_dir, "Input1_ImmunePanel_ALL_normalized_data_LN.xlsx", sep = "/")
group_xlsx <- paste(in_dir, "Input2_final_grouping_v1.xlsx", sep = "/")
result_xlsx <- paste(out_dir, "SD1_immune_ALL_normalized_log2_limma_LN_results.xlsx", sep = "/") ## Output: Supplementary Data 1
envol_plot <- paste(out_dir, "SF1_LN_High_vs_Low_EnhancedVolcanol_Immune_Only.png", sep = "/") ## Output: Supplementary Figure 1A

expr_df <- read_excel(expr_xlsx, sheet="log2")
expr_df <- as.data.frame(expr_df)
rownames(expr_df) <- expr_df$pid
expr_df$pid <- NULL
expr_df <- t(expr_df)
expr_pid <- colnames(expr_df)

group_df <- read_excel(group_xlsx)
group_df <- as.data.frame(group_df)
rownames(group_df) <- group_df$pid
group_pid <- rownames(group_df)
used_pid <- intersect(group_pid,expr_pid)


expr_df <- expr_df[,used_pid]
head(expr_df)
group_df <- group_df[used_pid,]
head(group_df)

gfac <- factor(group_df$group, levels=c("medium","low","high"))
gdes <- model.matrix(~gfac)

head(gdes)
dim(gdes)
con <- makeContrasts(test_diff = gfachigh-gfaclow, levels=gdes)
fit <- lmFit(expr_df, gdes)
fit2 <- contrasts.fit(fit, con)
fit2 <- eBayes(fit2)
top <- topTable(fit2, coef="test_diff", number=nrow(fit2), adjust.method="fdr")
print("test")
head(top)
write.xlsx(top, file=result_xlsx)
png(envol_plot, res=300, width=9, heigh=6, units='in')
EnhancedVolcano(top, lab=rownames(top), x="logFC", y="adj.P.Val", pCutoff=0.10, FCcutoff=1, xlim=c(-2,2), ylim=c(0, 1.5),
		title="TNBC tumor-negative lymph node (immune panel)", subtitle="Differential expression (High vs Low)", drawConnectors = FALSE,
		legendPosition="right",caption = "FC cutoff=1.0\t Adjusted p-value cutoff=0.10",
		hline = c(0.05, 0.01, 10e-3, 10e-5),)
gar <- dev.off()

