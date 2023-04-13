# This script is designed to merge and analyze the least-one correlation results
# Weihua Guo, Ph.D.,
# 10/22/2019

rm(list = ls())

suppressMessages(library(Hmisc))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(tibble))
suppressMessages(library(stringr))

rowDesStat = function(Df) {
  statDf = data.frame(matrix(ncol = 3, nrow = nrow(Df)))
  rownames(statDf) = rownames(Df)
  colnames(statDf) = c("pco005CTS", "pco001CTS", "Mean")
  statDf$pco005CTS = rowSums(Df<=0.05)
  statDf$pco001CTS = rowSums(Df<=0.01)
  statDf$Mean = rowMeans(Df)
  statDf$Median = apply(Df, 1, median)
  statDf$SD = apply(Df, 1, sd)
  statDf$SE = apply(Df, 1, sd)/sqrt(ncol(Df))
  statDf$max = apply(Df, 1, max)
  statDf$min = apply(Df, 1, min)
  return(statDf)
}

corMerge = function(groi, xLab, yLab, dataDir, csvSave = FALSE, rdsSave = TRUE) {
  st = Sys.time()
  filePattern = paste(groi, xLab, "vs", yLab, sep = "_")
  
  allFlatFiles = list.files(dataDir, pattern = "flat.csv")
  usedFiles = allFlatFiles[str_detect(allFlatFiles, filePattern)]
  
  mergeCor = NULL
  mergeP = NULL
  for (ifile in usedFiles) {
    message(ifile)
    tmpDf = read.csv(paste(dataDir, ifile, sep = ""), header = TRUE, row.names = 1)
    tmpClDf = tmpDf[tmpDf$row != tmpDf$column,]
    rownames(tmpClDf) = paste0(tmpClDf$row, "_vs_", tmpClDf$column)
    
    mergeCor = cbind(mergeCor, tmpClDf$cor)
    mergeP = cbind(mergeP, tmpClDf$p)
    
    # print(head(tmpClDf))
    
  }
  rownames(mergeCor) = rownames(tmpClDf)
  rownames(mergeP) = rownames(tmpClDf)
  colnames(mergeCor) = str_split_fixed(usedFiles, "_rem", n = 2)[,1]
  colnames(mergeP) = str_split_fixed(usedFiles, "_rem", n = 2)[,1]
  
  statP = rowDesStat(mergeP)
  mergeRes = list("r" = mergeCor, "P" = mergeP, "stat" = statP)
  
  if (csvSave) {
    write.csv(mergeCor, file = paste(dataDir, groi, "_merged_pcc_", xLab, "vs", yLab, ".csv", sep = ""))
    write.csv(mergeP, file = paste(dataDir, groi, "_merged_pcc_pvalue_", xLab, "vs", yLab, ".csv", sep = ""))
    write.csv(statP, file = paste(dataDir, groi, "_descriptive_stat_", xLab, "vs", yLab, ".csv", sep = ""))
  }
  
  if (rdsSave) {saveRDS(mergeRes, file = paste(dataDir, groi, "_merged_", xLab, "vs", yLab, ".RDS", sep = ""))}
  print(Sys.time()-st)
  return(mergeRes)
}

download_dir <- "..."# NOTE: Please change the ... into the directory you were saving the data

in_dir <- paste(download_dir, "FigS2/", sep = "/")
out_dir <- paste(download_dir, "FigS2/", sep = "/")

xLab = "LN"
yLab = "LN"

lowMres = corMerge("low", xLab, yLab, in_dir)
highMres = corMerge("high", xLab, yLab, in_dir)

xr = lowMres$r
yr = highMres$r
xrow = rownames(xr)
yrow = rownames(yr)
if (all(xrow == yrow)) {
  st = Sys.time()
  ttRes = data.frame(matrix(ncol = 6, nrow = length(xrow)))
  colnames(ttRes) = c("lowMean", "highMean", "tstat", "pval", "lCI", "hCI")
  rownames(ttRes) = xrow
  ict = 0
  for (irow in xrow) {
    tmpRes = t.test(lowMres$r[irow,], highMres$r[irow,])
    ttRes[irow, "tstat"] = tmpRes$statistic[[1]]
    ttRes[irow, "pval"] = tmpRes$p.value[[1]]
    ttRes[irow, "lowMean"] = tmpRes$estimate[[1]]
    ttRes[irow, "highMean"] = tmpRes$estimate[[2]]
    ttRes[irow, "lCI"] = tmpRes$conf.int[[1]]
    ttRes[irow, "hCI"] = tmpRes$conf.int[[2]]

    ict = ict + 1
    if (ict%%9000 == 0) {
      message(irow, ict)
    }
  }
  print(Sys.time()-st)
}

finalRes = ttRes
finalRes$MeanDiff = finalRes$highMean - finalRes$lowMean
finalRes$lowP005 = lowMres$stat[,"pco005CTS"]
finalRes$lowP001 = lowMres$stat[,"pco001CTS"]
finalRes$highP005 = highMres$stat[,"pco005CTS"]
finalRes$highP001 = highMres$stat[,"pco001CTS"]

write.csv(finalRes, file = paste(dataDir, "ttest_merge_results_", xLab, "vs", yLab, ".csv", sep = ""))
saveRDS(finalRes, file = paste(dataDir, "ttest_merge_results_", xLab, "vs", yLab, ".RDS", sep = ""))
