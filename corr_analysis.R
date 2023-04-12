library(dplyr)
library(readxl)
library(Hmisc)
library(xlsx)

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
out_dir <- paste(download_dir, "Fig2", sep = "/")

x_file <- paste(in_dir, "Input1_ImmunePanel_ALL_normalized_data_LN.xlsx", sep = "/")
x_lab <- "LN"
group_file <- paste(in_dir, "Input2_final_grouping_v1.xlsx", sep = "")

################################
# Parameter section
igoi <- "low" # Immune Group Of Interests ## NOTE: Need to run for low and high separatively.
corr_method <- "pearson"
data_type <- "immune_only"
res_type <- "xself" # avsb or aself or bself

#################
# save results to tables
table_filenames <- c("r","p","fm") # fm: flat_matrix
table_save <- vector(mode="list", length=length(table_filenames))
names(table_save) <- table_filenames
table_save["r"] <- paste(out_dir,"corr_coef",igoi,corr_method,data_type,res_type,".csv",sep="/_") # Input for heatmap_plotter.R
table_save["p"] <- paste(out_dir,"p_value",igoi,corr_method,data_type,res_type,".csv",sep="/_") # Input for heatmap_plotter.R
table_save["fm"] <- paste(out_dir,"flat_corr_res",igoi,corr_method,data_type,res_type,".csv",sep="/_") # Output: Supplementary Data 3

# Read x matrix: LN
x_df <- read_excel(x_file)
x_df <- as.data.frame(x_df)
rownames(x_df) <- x_df$pid
x_df$pid <- NULL
x_df <- t(x_df)
# head(x_df)

#Read grouping info.
group_df <- read_excel(group_file)
group_df <- as.data.frame(group_df)
rownames(group_df) <- group_df$pid
# head(group_df)
# Filter by group
igoi_mask <- group_df$group == igoi
# igoi_mask
igpid <- rownames(group_df[igoi_mask,])
# length(igpid)

x_igoi <- x_df[,igpid]
# dim(x_igoi)
# Add substring to gene names for labeling data resource
rownames(x_igoi) <- paste(rownames(x_igoi),"_",x_lab,sep="")
# head(x_igoi)

cor_res <- rcorr(t(x_igoi), type=corr_method) ##!!! THIS DID ALL THE POSSIBLE correlations

xrow <- length(rownames(x_igoi))

res_r <- cor_res$r[1:xrow,1:xrow]
res_p <- cor_res$P[1:xrow,1:xrow]
# Remove the upper triangle for self correlated correlations
res_rl <- res_r
res_pl <- res_p
res_rl[lower.tri(res_rl)] <- NA
res_pl[lower.tri(res_rl)] <- NA
res_flat <- flat_cor_mat(res_rl, res_pl)
res_flat <- res_flat[complete.cases(res_flat[,3]),]

cat("Table save files:\n")
print(table_save)
write.csv(res_r, file=table_save[["r"]])
write.csv(res_p, file=table_save[["p"]])
write.csv(res_flat, file=table_save[["fm"]])
