# This script is for visualization of the least one correlation 
# Weihua Guo, Ph.D.
# 10/26/2019

rm(list = ls())

suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(readxl))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(stringr))
suppressMessages(library(EnhancedVolcano))

download_dir <- "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir <- paste(download_dir, "FigS2/", sep = "/")
out_dir <- paste(download_dir, "FigS2/", sep = "/")

lowRDS = "low_merged_LNvsLN.RDS"
highRDS = "high_merged_LNvsLN.RDS"

ttRDS = "ttest_merge_results_LNvsLN.RDS"
valOI = c("Mean", "pco005CTS", "pco001CTS", "Median", "max", "min")

st = Sys.time()
message("Start to read results...")
lowRes = readRDS(paste(in_dir, lowRDS, sep = ""))
highRes = readRDS(paste(in_dir, highRDS, sep = ""))
ttRes = readRDS(paste(in_dir, ttRDS, sep = ""))
print(Sys.time() - st)

pdf(paste(out_dir, "pcc_pval_volcano.pdf", sep = ""), width = 19.2, heigh = 10.8) # Output: Fig S2
EnhancedVolcano(ttRes,
		lab = rownames(ttRes),
		x = "MeanDiff",
		y = "pval",
		FCcutoff = 1.5,
		title = NULL,
		subtitle = NULL,
		xlab = "Delta PCC (High-Low)",
		legendPosition = "right")
gar = dev.off()
