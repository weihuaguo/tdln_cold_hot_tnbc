rm(list = ls())
suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
data_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/ihc/MC_IL4_mDC_Th2_Panel/roi_analysis"
qp_annot_file <- "20240121_2121_20ROIs_annotations.tsv"
qp_annot_file <- "20240123_2250_30ROIs_annotations.tsv"
smp_annot_file <- "sample_collection_231128.xlsx"

annot_df <- read.table(paste(data_dir, qp_annot_file, sep = "/"), check.names = F, header = T, sep = "\t")
annot_df[is.na(annot_df)] <- 0
annot_df$sid <- str_split_fixed(annot_df$Image, "_b", n = 2)[,1]
print(head(annot_df))

smp_df <- as.data.frame(read_excel(paste(data_dir, smp_annot_file, sep = "/")))
smp_df$sid <- str_replace_all(smp_df$`Slide ID`, "_", "-")
smp_df <- smp_df[smp_df$sid != "NA",]
print(head(smp_df))

merge_df <- merge(annot_df, smp_df, by = "sid", all.x = T)
#write.csv(merge_df, paste (data_dir, "20240121_2121_20ROI_merge_annotation.csv", sep = "/"))
write.csv(merge_df, paste (data_dir, "20240123_2250_30ROI_merge_annotation.csv", sep = "/"))
print(head(merge_df))

sum_df <- merge_df %>%
	group_by(sid, Cohort) %>%
	summarize(n_roi = n()/3,
		  area = sum(`Area µm^2`),
		  `IL4+CD117+` =  sum(`Num IL4: CD117`),
		  `CD117+IL4+` =  sum(`Num CD117: IL4`),
		  `IL4+` =  sum(`Num IL4`),
		  `CD117+` =  sum(`Num CD117`),
	)
print(sum_df)
write.csv(sum_df, paste (data_dir, "20240123_2250_30ROI_merge_annotation_slide_wise.csv", sep = "/"))

