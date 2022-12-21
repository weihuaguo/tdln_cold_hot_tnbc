# AUC for cell densities in TDLN of cold vs hot
# Weihua Guo, Ph.D.
# 22.12.21

rm(list = ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(precrec))
suppressMessages(library(readxl))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))

data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/ihc"

input_df <- as.data.frame(read_excel(paste(data_dir, "for_auc_input.xlsx", sep = "/")))
print(input_df)
print(colnames(input_df))

cat("For density multiple tests...\n")
use_df <- input_df[,c("Cohort", colnames(input_df)[7:ncol(input_df)])]
use_df$label <- ifelse(use_df$Cohort == "COLD", 1, 0)
print(use_df)


scores <- join_scores(use_df[,2:(ncol(use_df)-1)])
msmdat <- mmdata(scores, use_df$label, modnames = colnames(use_df)[2:(ncol(use_df)-1)])
msccurves <- evalmod(msmdat)
auc_curves <- auc(msccurves)
write.csv(auc_curves, paste(data_dir, "density_auc_results.csv", sep = "/"))

gg <- autoplot(msccurves, curvetype = c("ROC"))
ggsave(paste(data_dir, "density_auc_roc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)

gg <- autoplot(msccurves, curvetype = c("PRC"))
ggsave(paste(data_dir, "density_auc_prc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)

cat("For ratio multiple tests...\n")
use_df <- input_df[,c("Cohort", colnames(input_df)[3:6])]
use_df$label <- ifelse(use_df$Cohort == "COLD", 1, 0)
print(use_df)
scores <- join_scores(use_df[,2:(ncol(use_df)-1)])
msmdat <- mmdata(scores, use_df$label, modnames = colnames(use_df)[2:(ncol(use_df)-1)])
msccurves <- evalmod(msmdat)
auc_curves <- auc(msccurves)
write.csv(auc_curves, paste(data_dir, "ratio_auc_results.csv", sep = "/"))

gg <- autoplot(msccurves, curvetype = c("ROC"))
ggsave(paste(data_dir, "ratio_auc_roc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)

gg <- autoplot(msccurves, curvetype = c("PRC"))
ggsave(paste(data_dir, "ratio_auc_prc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)

