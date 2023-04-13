# This script is designed to generate least-one correlation coeffients
# Weihua Guo, Ph.D.,
# 10/22/2019

rm(list = ls())

suppressMessages(library(Hmisc))
suppressMessages(library(dplyr))
suppressMessages(library(readxl))
suppressMessages(library(tibble))

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

download_dir <- "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir <- paste(download_dir, "input_file", sep = "/")
out_dir <- paste(download_dir, "FigS2", sep = "/")

dataDir = "/home/weihua/mnts/group_plee/Weihua/nanostring_results/outline_and_lists/"
resDir = "/home/weihua/mnts/group_plee/Weihua/nanostring_results/least_one_results/"
xFile = "Input1_ImmunePanel_ALL_normalized_data_LN.xlsx"
xLabs = "LN"
yFile = "Input1_ImmunePanel_ALL_normalized_data_LN.xlsx"
yLabs = "LN"

groupXlsx = "Input2_final_grouping_v1.xlsx"

###########
# Frequent changes
groi = "high" # NOTE: Need to run for "low" and "high" separatively.

xDf = read_excel(paste(in_dir, xFile, sep = ""))
xDf = column_to_rownames(xDf, "pid")
colnames(xDf) = paste0(colnames(xDf), "_", xLabs)
yDf = read_excel(paste(in_dir, yFile, sep = ""))
yDf = column_to_rownames(yDf, "pid")
colnames(yDf) = paste0(colnames(yDf), "_", yLabs)

groupDf = read_excel(paste(in_dir, groupXlsx, sep = ""))

pidoi = groupDf[groupDf$group == groi, "pid"]
xoiDf = xDf[pidoi[["pid"]],]
yoiDf = yDf[pidoi[["pid"]],]

if (sum(pidoi[["pid"]] %in% rownames(xDf))[1] != length(pidoi[["pid"]])) {
  stop("MISSING PID!!!")
} else {message("Patient of interests are extracted.")}

for (ipid in pidoi[["pid"]]) {
  st = Sys.time()
  message("Start to run the correlation analysis!")
  message("\t",ipid)
  subXoiDf = xoiDf[-which(rownames(xoiDf) == ipid),]
  subYoiDf = yoiDf[-which(rownames(xoiDf) == ipid),]

  cor_res = rcorr(x = as.matrix(subXoiDf), y = as.matrix(subYoiDf))
  if (xLabs == yLabs) {
  corR = cor_res$r[1:ncol(xDf), 1:ncol(yDf)]
  corP = cor_res$P[1:ncol(xDf), 1:ncol(yDf)]
  } else {
    corR = cor_res$r[1:ncol(xDf), (ncol(yDf)+1):ncol(cor_res$r)]
    corP = cor_res$P[1:ncol(xDf), (ncol(yDf)+1):ncol(cor_res$P)]
  }
  flatRes = flat_cor_mat(corR, corP)
  
  message("\tStart to output results...")
  rPrefix = paste(resDir, groi, "_", xLabs, "_vs_", yLabs, "_", ipid, "_removed_r.csv", sep = "")
  pPrefix = paste(resDir, groi, "_", xLabs, "_vs_", yLabs, "_", ipid, "_removed_p.csv", sep = "")
  fPrefix = paste(resDir, groi, "_", xLabs, "_vs_", yLabs, "_", ipid, "_removed_flat.csv", sep = "")
  
  write.csv(corR, file = rPrefix)
  write.csv(corP, file = pPrefix)
  write.csv(flatRes, file = fPrefix)
  print(Sys.time() - st)
}
