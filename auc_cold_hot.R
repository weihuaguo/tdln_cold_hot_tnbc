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

dis_mast <- c(123.3438654, 639.105862, 321.2550346, 78.50367433, 369.9601241, 310.4804941, 204.7372405, 31.80297459, 
	      70.14630692, 42.56299425, 132.7860555, 52.88749536, 41.9307874, 60.78878487)
val_mast <- c(50.00714985, 43.16583183, 39.03962527, 48.64217744, 282.3807063, 27.34932313, 159.2904235, 113.4510706, 
	      152.4983617, 164.4171712, 443.8816994, 58.74099658, 141.517457, 138.1741809, 110.0876861, 352.1657489, 431.9875945, 225.6056757, 356.7502554, 277.8699265)
dis_lab <- c(1,1,1,1,1,1,1,0,0,0,0,0,0,0)
val_lab <- c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)
all_mast <- c(dis_mast, val_mast)
all_lab <- c(dis_lab, val_lab)

#use_score <- join_scores(dis_mast, val_mast, all_mast)
print("FuckR")

#use_label <- join_labels(dis_lab, val_lab, all_lab)

msmdat <- mmdata(dis_mast, dis_lab, modnames = c("Discovery"))
msccurves <- evalmod(msmdat)
dis_auc_curves <- auc(msccurves)
write.csv(dis_auc_curves, paste(data_dir, "dis_mast_auc_results.csv", sep = "/"))

gg <- autoplot(msccurves, curvetype = c("ROC"))
ggsave(paste(data_dir, "dis_mast_auc_roc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)
pgg <- ggplot_build(gg)
dis_roc <- pgg$data[[1]]
dis_roc$cohort <- "Discovery"

gg <- autoplot(msccurves, curvetype = c("PRC"))
ggsave(paste(data_dir, "dis_mast_auc_prc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)

msmdat <- mmdata(val_mast, val_lab, modnames = c("Validation"))
msccurves <- evalmod(msmdat)
val_auc_curves <- auc(msccurves)
write.csv(val_auc_curves, paste(data_dir, "val_mast_auc_results.csv", sep = "/"))

gg <- autoplot(msccurves, curvetype = c("ROC"))
ggsave(paste(data_dir, "val_mast_auc_roc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)
pgg <- ggplot_build(gg)
val_roc <- pgg$data[[1]]
val_roc$cohort <- "Validation"

gg <- autoplot(msccurves, curvetype = c("PRC"))
ggsave(paste(data_dir, "val_mast_auc_prc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)


msmdat <- mmdata(all_mast, all_lab, modnames = c("Validation"))
msccurves <- evalmod(msmdat)
all_auc_curves <- auc(msccurves)
write.csv(all_auc_curves, paste(data_dir, "all_mast_auc_results.csv", sep = "/"))

gg <- autoplot(msccurves, curvetype = c("ROC"))
ggsave(paste(data_dir, "all_mast_auc_roc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)
pgg <- ggplot_build(gg)
all_roc <- pgg$data[[1]]
all_roc$cohort <- "Both"

gg <- autoplot(msccurves, curvetype = c("PRC"))
ggsave(paste(data_dir, "all_mast_auc_prc.png", sep = "/"), dpi = 300, width = 5, height = 4.2)

merge_roc <- rbind(dis_roc, val_roc)
merge_roc <- rbind(merge_roc, all_roc)

merge_gg <- ggplot(merge_roc, aes(x = x, y = y, color = cohort, group = cohort)) +
	geom_line() +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
	labs(x = "1 - Specificity", y = "Sensitivity", color = "Dataset") +
	theme_bw()
ggsave(paste(data_dir, "merged_mast_auc_roc.png", sep = "/"), dpi = 300, width = 5.4, height = 4.2)

q(save = "no")
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

