rm(list = ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))

download_dir <- "/home/weihua/mnts/smb_plee/Group/weihua/tdln_cold_hot_data_hub/ihc"# NOTE: Please change the ... into the directory you were saving the data

panel_dir <- paste(download_dir, "mDC_Th_Panel", sep = "/")
fr_dir <- paste(panel_dir, "fixed_radius", sep = "/")
nn_dir <- paste(panel_dir, "nearest_neighbors", sep = "/")

frs <- c(10,20,30,40,50)
smp_df <- as.data.frame(read_excel(paste(panel_dir, "sample_annotation_ihc_discovery.xlsx", sep = "/")))

fr_flag <- T
merge_fr_flag <- F
vis_fr_flag <- T
nn_flag <- F

if (fr_flag) {
	if (merge_fr_flag) {
		for (ifr in frs) {
			cat("\t", ifr, "\n")
			all_res <- list.files(fr_dir, pattern = paste("fr", ifr, sep = ""))
			print(all_res)
			c <- 0
			for (ifile in all_res) {
				cat("\t\t", ifile, "\n")
				df <- read.csv(paste(fr_dir, ifile, sep = "/"), row.names = 1, check.names = F)
	#			print(head(df))
				mdc_df <- df[df$phenotype == "mDC",]
				mdc_df$`r_CD8+` <- mdc_df$`CD8+`/sum(df$phenotype == "CD8+")
				mdc_df$`r_Th2` <- mdc_df$Th2/sum(df$phenotype == "Th2")
				mdc_df$`r_Th1` <- mdc_df$Th1/sum(df$phenotype == "Th1")

				mdc_df$th21_ratio <- (mdc_df$Th2+0.1)/(mdc_df$Th1+0.1)
				mdc_df$th21_diff <- mdc_df$Th2-mdc_df$Th1

				mdc_df$th2_larger_flag <- ifelse(mdc_df$Th2 > mdc_df$Th1, "Th2L", "Th1LE")
				mdc_df$both0_flag <- ifelse(mdc_df$Th2==0 & mdc_df$Th1==0, "Both0", "N")
				mdc_df$th1_0_flag <- ifelse(mdc_df$Th1==0 & mdc_df$Th2>0, "Only_Th1_0", "N")

				mdc_th_df <- mdc_df[,c("Th1", "Th2")]
				mdc_df$th_max_flag <- colnames(mdc_th_df)[apply(mdc_th_df, 1, which.max)]
				mdc_df$th_flag <- ifelse(mdc_df$both0_flag == "Both0", "No Th cell around",
							 ifelse(mdc_df$Th2 == mdc_df$Th1, "Th2 = Th1", mdc_df$th_max_flag)
				)
				mdc_df$total_th2 <- sum(df$phenotype == "Th2")
				mdc_df$total_th1 <- sum(df$phenotype == "Th1")

# TODO: based on relative Th2 and Th1 numbers divided by total Th2/Th1 cell numbers

	#			mdc_t_df <- mdc_df[,c("CD8+", "Th1", "Th2")]
	#			mdc_df$t_flag <- colnames(mdc_t_df)[apply(mdc_t_df, 1, which.max)]
	#			mdc_df$abs_def <- ifelse(mdc_df$both0_flag == "Both0", "No T cells", mdc_df$t_flag)

	#			mdc_t_df <- mdc_df[,c("r_CD8+", "r_Th1", "r_Th2")]
	#			mdc_df$r_t_flag <- colnames(mdc_t_df)[apply(mdc_t_df, 1, which.max)]
	#			mdc_df$rel_def <- ifelse(mdc_df$both0_flag == "Both0", "No T cells", mdc_df$r_t_flag)
			
	#			print(head(mdc_df))

				if (c == 0) {
					merge_mdc_df <- mdc_df
				} else {
					merge_mdc_df <- rbind(merge_mdc_df, mdc_df)
				}
				c <- c+1
			}
			write.csv(merge_mdc_df, paste(fr_dir, "/merge_fr_", ifr, "_cellwise.csv", sep = ""))
		}
	}

	if (vis_fr_flag) {
		for (ifr in frs) {
			cat(ifr, "\n")
			fr_df <- read.csv(paste(fr_dir, "/merge_fr_", ifr, "_cellwise.csv", sep = ""), check.names = F, row.names = 1)
			print(head(fr_df))
			fr_df <- fr_df %>% group_by(Image) %>% mutate(total_mdc = n())

			sum_df <- fr_df %>% 
				group_by(Image) %>%
				summarize(Th1 = sum(Th1), Th2 = sum(Th2), total_th1 = mean(total_th1), total_th2 = mean(total_th2))
			sum_df$th21_ratio <- sum_df$Th2/sum_df$Th1
			sum_df$sld_id <- str_replace_all(str_split_fixed(sum_df$Image, "x", n = 2)[,1], "-", "_")
			sum_df <- merge(sum_df, smp_df, by = "sld_id", all.x = T)
			sum_df$r_th1 <- sum_df$Th1/sum_df$total_th1
			sum_df$r_th2 <- sum_df$Th2/sum_df$total_th2
			sum_df$r_th21_ratio <- sum_df$r_th2/sum_df$r_th1
#			print(head(sum_df))
#			q(save = "no")

			box_gg <- ggplot(sum_df, aes(x = cohort, y = th21_ratio, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Th2/Th1 ratio") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_sum_th21_ratio.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = r_th21_ratio, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Relative Th2/Th1 ratio") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_sum_rel_th21_ratio.png", sep = ""), dpi = 300, width = 3, height = 6)


			cat("\tAll mDCs\n")
			mdc_cts_df <- fr_df %>%
				group_by(Image, th_flag) %>%
				summarize(n_mdc = n(), total_mdc = mean(total_mdc))
			mdc_cts_df$sld_id <- str_replace_all(str_split_fixed(mdc_cts_df$Image, "x", n = 2)[,1], "-", "_")
			mdc_cts_df <- as.data.frame(mdc_cts_df)
			mdc_cts_df$rel_mdc <- mdc_cts_df$n_mdc/mdc_cts_df$total_mdc*100
#			print(head(mdc_cts_df))
#			print(head(smp_df))

			spr_df <- spread(mdc_cts_df[,c("th_flag", "rel_mdc", "sld_id")], "th_flag", "rel_mdc")
			spr_df[is.na(spr_df)] <- 0.0
			gath_df <- gather(spr_df, "th_flag", "rel_mdc", unique(mdc_cts_df$th_flag))
			gath_df <- merge(gath_df, smp_df, by = "sld_id", all.x = T)
#			print(head(gath_df))

			box_gg <- ggplot(gath_df, aes(x = cohort, y = rel_mdc, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				facet_wrap(~th_flag, nrow = 1, scales = "free") +
				labs(title = paste("Fixed radius:", ifr), y = "Relative to all mDCs") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_total_mdc_percentage.png", sep = ""), dpi = 300, width = 6, height = 6)

			cat("\tmDCs excluding the no Th cell around\n")
			mdc_cts_df <- fr_df %>%
				group_by(Image, th_flag) %>%
				summarize(n_mdc = n(), total_mdc = mean(total_mdc))
			mdc_cts_df$sld_id <- str_replace_all(str_split_fixed(mdc_cts_df$Image, "x", n = 2)[,1], "-", "_")
			mdc_cts_df <- as.data.frame(mdc_cts_df)
#			print(head(mdc_cts_df))
#			print(head(smp_df))

			spr_df <- spread(mdc_cts_df[,c("th_flag", "n_mdc", "sld_id", "total_mdc")], "th_flag", "n_mdc")
			spr_df[is.na(spr_df)] <- 0.0
			spr_df$total_mdc <- spr_df$total_mdc-spr_df$`No Th cell around`
			print(head(spr_df))
			gath_df <- gather(spr_df, "th_flag", "n_mdc", unique(mdc_cts_df$th_flag))
			gath_df$rel_mdc <- gath_df$n_mdc/gath_df$total_mdc*100
			gath_df <- merge(gath_df, smp_df, by = "sld_id", all.x = T)
#			print(head(gath_df))

			box_gg <- ggplot(gath_df, aes(x = cohort, y = rel_mdc, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				facet_wrap(~th_flag, nrow = 1, scales = "free") +
				labs(title = paste("Fixed radius:", ifr), y = "Realtive to mDCs with Th around") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_tharound_mdc_percentage.png", sep = ""), dpi = 300, width = 6, height = 6)

		}
	}
}

