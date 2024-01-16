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
frs <- c(10,20,30,40,50,75,100,125,150,200)
frs <- c(10,20,30,40,50,75,100)


smp_df <- as.data.frame(read_excel(paste(panel_dir, "sample_annotation_ihc_discovery.xlsx", sep = "/")))

fr_flag <- T
merge_fr_flag <- F
vis_fr_flag <- T
nn_flag <- F
merge_nn_flag <- F
vis_nn_flag <- F
th21_nn_flag <- F
th21_vis_nn_flag <- F

if (th21_nn_flag) {
	if (merge_nn_flag) {
		c <- 0
		all_res <- list.files(nn_dir, pattern = "nearest_neighbor_distances")
		for (ifile in all_res) {
			cat(ifile, "\n")
			df <- read.csv(paste(nn_dir, ifile, sep = "/"), row.names = 1, check.names = F)
			th21_df <- df[df$phenotype %in% c("Th2", "Th1", "CD8+"),]
			print(head(th21_df))
			mdc_df <- as.data.frame(matrix(ncol = length(frs), nrow = nrow(th21_df)))
			rownames(mdc_df) <- rownames(th21_df)
			colnames(mdc_df) <- paste("mdc_within_", frs, sep = "")
			for (ifr in frs) {
				mdc_df[,paste("mdc_within_", ifr, sep = "")] <- ifelse(th21_df$mDC <= ifr, paste("Within", ifr, "um"), "No")
			}

			mdc_df <- merge(th21_df, mdc_df, by = "row.names", all.x = T)
			print(head(mdc_df))

			if (c == 0) {
				merge_mdc_df <- mdc_df
			} else {
				merge_mdc_df <- rbind(merge_mdc_df, mdc_df)
			}
			c <- c+1
		}
		write.csv(merge_mdc_df, paste(nn_dir, "/merge_nn_th21_fr2mdc_df_cellwise.csv", sep = ""))
	}
	if (th21_vis_nn_flag) {
		mdc_df <- read.csv(paste(nn_dir, "/merge_nn_th21_fr2mdc_df_cellwise.csv", sep = ""), row.names = 1, check.names = F)
		mdc_df <- mdc_df %>% group_by(Image, phenotype) %>% mutate(n_th = n())
		gath_df <- gather(mdc_df, "dist_cutoff", "status", colnames(mdc_df)[str_detect(colnames(mdc_df), "_within_")])
		print(head(gath_df))
		sum_df <- gath_df %>%
			group_by(Image, phenotype, dist_cutoff, status) %>%
			summarize(n = n(), n_th = mean(n_th))
		sum_df$rel <- sum_df$n/sum_df$n_th*100
		sum_df$sld_id <- str_replace_all(str_split_fixed(sum_df$Image, "x", n = 2)[,1], "-", "_")
		sum_df <- merge(sum_df, smp_df, by = "sld_id", all.x = T)
		sum_df <- sum_df[sum_df$status != "No",]
		sum_df$dist_cutoff <- factor(sum_df$dist_cutoff, levels = paste("mdc_within_", frs, sep = ""))
		write.csv(sum_df, paste(nn_dir, "/merge_nn_th21_fr2med_df_summary.csv", sep = ""))

		box_gg <- ggplot(sum_df, aes(x = cohort, y = rel, color = cohort)) +
			geom_boxplot() +
			geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			facet_grid(phenotype~dist_cutoff, scales = 'free') +
			labs(y = "Th percentage with nearest mDC in x um") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_th_pct_with_mdc_in_fr.png", sep = ""), dpi = 300, width = 16, height = 12)

	}
}

if (nn_flag) {
	if (merge_nn_flag) {
		c <- 0
		all_res <- list.files(nn_dir, pattern = "nearest_neighbor_distances")
		for (ifile in all_res) {
			cat(ifile, "\n")
			df <- read.csv(paste(nn_dir, ifile, sep = "/"), row.names = 1, check.names = F)
			mdc_df <- df[df$phenotype == "mDC",]
			mdc_df$th21_flag <- ifelse(mdc_df$Th2 > mdc_df$Th1, "Th1", "Th2")
			tdist_df <- mdc_df[,c("CD8+", "Th2", "Th1")]
			adist_df <- mdc_df[,c("CD8+", "Th2", "Th1", "mDC")]

			mdc_df$t_flag <- colnames(tdist_df)[apply(tdist_df, 1, which.min)]
			mdc_df$a_flag <- colnames(adist_df)[apply(adist_df, 1, which.min)]
			if (c == 0) {
				merge_mdc_df <- mdc_df
			} else {
				merge_mdc_df <- rbind(merge_mdc_df, mdc_df)
			}
			c <- c+1
		}
		write.csv(merge_mdc_df, paste(nn_dir, "/merge_nn_mdc_df_cellwise.csv", sep = ""))
	}

	if (vis_nn_flag) {
		mdc_df <- read.csv(paste(nn_dir, "/merge_nn_mdc_df_cellwise.csv", sep = ""), row.names = 1, check.names = F)
		print(head(mdc_df))

		dist_ratio_df <- mdc_df
		dist_ratio_df$th21_dr <- dist_ratio_df$Th2/dist_ratio_df$Th1

		dist_ratio_df$sld_id <- str_replace_all(str_split_fixed(dist_ratio_df$Image, "x", n = 2)[,1], "-", "_")
		gath_df <- merge(dist_ratio_df, smp_df, by = "sld_id", all.x = T)

		smpl_avg_df <- gath_df %>%
			group_by(cohort, sld_id) %>%
			summarise(n = n(),
				  mean = mean(th21_dr, na.rm = T),
				  median = median(th21_dr, na.rm = T),
				  sd = sd(th21_dr, na.rm = T),
				  min = min(th21_dr, na.rm = T),
				  max = max(th21_dr, na.rm = T),
				  se = sd/sqrt(n),
				  mean_se_upper = mean+se,
				  mean_se_lower = mean-se)
		write.csv(smpl_avg_df, paste(nn_dir, "/merge_nn_descriptive_stat_samplewise_th21_dist_ratio.csv", sep = ""))

		box_gg <- ggplot(smpl_avg_df, aes(x = cohort, y = mean, color = cohort)) +
			geom_boxplot(outlier.shape = NA) +
			geom_point(position = position_jitterdodge(dodge.width = 0.7), alpha = 0.666666) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			labs(y = "Distance ratio between nearest Th2 and Th1, sample wise") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_th21_dist_ratio_samplewise.png", sep = ""), dpi = 300, width = 3, height = 6)

		cell_avg_df <- gath_df %>%
			group_by(cohort) %>%
			summarise(n = n(),
				  mean = mean(th21_dr, na.rm = T),
				  median = median(th21_dr, na.rm = T),
				  sd = sd(th21_dr, na.rm = T),
				  min = min(th21_dr, na.rm = T),
				  max = max(th21_dr, na.rm = T),
				  se = sd/sqrt(n),
				  mean_se_upper = mean+se,
				  mean_se_lower = mean-se)
		write.csv(cell_avg_df, paste(nn_dir, "/merge_nn_descriptive_stat_cellwise_th21_dist_ratio.csv", sep = ""))


		box_gg <- ggplot(gath_df, aes(x = cohort, y = th21_dr, color = cohort)) +
			geom_boxplot(outlier.shape = NA) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			labs(y = "Distance ratio between nearest Th2 and Th1, cellwise") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_th21_dist_ratio_cellwise.png", sep = ""), dpi = 300, width = 3, height = 6)

		gath_df <- gather(mdc_df, "nearest_cell", "dist", c("mDC", "CD8+", "Th2", "Th1"))
		print(head(gath_df))
		gath_df$sld_id <- str_replace_all(str_split_fixed(gath_df$Image, "x", n = 2)[,1], "-", "_")
		gath_df <- merge(gath_df, smp_df, by = "sld_id", all.x = T)

		smpl_avg_df <- gath_df %>%
			group_by(cohort, nearest_cell, sld_id) %>%
			summarise(n = n(),
				  mean = mean(dist, na.rm = T),
				  median = median(dist, na.rm = T),
				  sd = sd(dist, na.rm = T),
				  min = min(dist, na.rm = T),
				  max = max(dist, na.rm = T),
				  se = sd/sqrt(n),
				  mean_se_upper = mean+se,
				  mean_se_lower = mean-se)
		write.csv(smpl_avg_df, paste(nn_dir, "/merge_nn_descriptive_stat_samplewise.csv", sep = ""))

		box_gg <- ggplot(smpl_avg_df, aes(x = cohort, y = mean, color = cohort)) +
			geom_boxplot(outlier.shape = NA) +
			geom_point(position = position_jitterdodge(dodge.width = 0.7), alpha = 0.666666) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			facet_wrap(~nearest_cell, nrow = 1, scales = "free") +
			labs(y = "Distance of nearest cells, sample_wise") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_dist_samplewise.png", sep = ""), dpi = 300, width = 6, height = 6)


		cell_avg_df <- gath_df %>%
			group_by(cohort, nearest_cell) %>%
			summarise(n = n(),
				  mean = mean(dist, na.rm = T),
				  median = median(dist, na.rm = T),
				  sd = sd(dist, na.rm = T),
				  min = min(dist, na.rm = T),
				  max = max(dist, na.rm = T),
				  se = sd/sqrt(n),
				  mean_se_upper = mean+se,
				  mean_se_lower = mean-se)
		write.csv(cell_avg_df, paste(nn_dir, "/merge_nn_descriptive_stat_cellwise.csv", sep = ""))

		box_gg <- ggplot(gath_df, aes(x = cohort, y = dist, color = cohort)) +
			geom_boxplot(outlier.shape = NA) +
#			geom_point(position = position_jitterdodge(dodge.width = 0.7), alpha = 0.001) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			facet_wrap(~nearest_cell, nrow = 1, scales = "free") +
			labs(y = "Distance of nearest cells, cell-wise") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_dist_cellwise.png", sep = ""), dpi = 300, width = 6, height = 6)

		a_sum_df <- mdc_df %>%
			group_by(Image) %>%
			mutate(n_mdc = n())
		a_sum_df <- a_sum_df %>%
			group_by(Image, a_flag) %>%
			summarize(n_a = n(), n_mdc = mean(n_mdc)) %>%
			mutate(rel_a = n_a/n_mdc*100)
		print(a_sum_df)
		a_sum_df$sld_id <- str_replace_all(str_split_fixed(a_sum_df$Image, "x", n = 2)[,1], "-", "_")
		sum_df <- merge(a_sum_df, smp_df, by = "sld_id", all.x = T)
		sum_df[is.na(sum_df)] <- 0.0
		box_gg <- ggplot(sum_df, aes(x = cohort, y = rel_a, color = cohort)) +
			geom_boxplot() +
			geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			facet_wrap(~a_flag, nrow = 1, scales = "free") +
			labs(y = "mDC percentage (%) (Th2, Th1, CD8+ and mDC)") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_sum_rel_mdc_all.png", sep = ""), dpi = 300, width = 6, height = 6)


		t_sum_df <- mdc_df %>%
			group_by(Image) %>%
			mutate(n_mdc = n())
		t_sum_df <- t_sum_df %>%
			group_by(Image, t_flag) %>%
			summarize(n_t = n(), n_mdc = mean(n_mdc)) %>%
			mutate(rel_t = n_t/n_mdc*100)
		print(t_sum_df)
		t_sum_df$sld_id <- str_replace_all(str_split_fixed(t_sum_df$Image, "x", n = 2)[,1], "-", "_")
		sum_df <- merge(t_sum_df, smp_df, by = "sld_id", all.x = T)
		sum_df[is.na(sum_df)] <- 0.0
		box_gg <- ggplot(sum_df, aes(x = cohort, y = rel_t, color = cohort)) +
			geom_boxplot() +
			geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			facet_wrap(~t_flag, nrow = 1, scales = "free") +
			labs(y = "mDC percentage (%) (Th2, Th1, CD8+ only)") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_sum_rel_mdc_t_only.png", sep = ""), dpi = 300, width = 4, height = 6)

		th21_sum_df <- mdc_df %>%
			group_by(Image) %>%
			mutate(n_mdc = n())
		th21_sum_df <- th21_sum_df %>%
			group_by(Image, th21_flag) %>%
			summarize(n_th21 = n(), n_mdc = mean(n_mdc)) %>%
			mutate(rel_th21 = n_th21/n_mdc*100)
		print(th21_sum_df)
		th21_sum_df$sld_id <- str_replace_all(str_split_fixed(th21_sum_df$Image, "x", n = 2)[,1], "-", "_")
		sum_df <- merge(th21_sum_df, smp_df, by = "sld_id", all.x = T)
		sum_df[is.na(sum_df)] <- 0.0
		box_gg <- ggplot(sum_df, aes(x = cohort, y = rel_th21, color = cohort)) +
			geom_boxplot() +
			geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
			stat_compare_means(aes(group = cohort), label = "p.format") +
			facet_wrap(~th21_flag, nrow = 1, scales = "free") +
			labs(y = "mDC percentage (%) (Th2, Th1 only)") +
			theme_bw()
		ggsave(paste(nn_dir, "/merge_nn_boxplot_sum_rel_mdc_th21_only.png", sep = ""), dpi = 300, width = 4, height = 6)


	}
}

if (fr_flag) {
	if (merge_fr_flag) {
		for (ifr in frs) {
			cat("\t", ifr, "\n")
			all_res <- list.files(fr_dir, pattern = paste("fr", ifr, "_", sep = ""))
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
				mdc_df$total_cd8t <- sum(df$phenotype == "CD8+")

				t_df <- df[df$phenotype %in% c("Th1", "Th2", "CD8+"),]
				t_df <- t_df[t_df$mDC>0,]
				mdc_df$around_th2 <- sum(t_df$phenotype == "Th2")
				mdc_df$around_th1 <- sum(t_df$phenotype == "Th1")
				mdc_df$around_cd8t <- sum(t_df$phenotype == "CD8+")
				
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
				summarize(sum_around_th1 = sum(Th1), sum_around_th2 = sum(Th2), sum_around_cd8t = sum(`CD8+`),
					  avg_around_th1 = mean(Th1), avg_around_th2 = mean(Th2), avg_around_cd8t = mean(`CD8+`),
					  true_around_th1 = mean(around_th1), true_around_th2 = mean(around_th2), true_around_cd8t = mean(around_cd8t),
					  total_th1 = mean(total_th1), total_th2 = mean(total_th2), total_cd8t = mean(total_cd8t))
			sum_df$sld_id <- str_replace_all(str_split_fixed(sum_df$Image, "x", n = 2)[,1], "-", "_")
			sum_df <- merge(sum_df, smp_df, by = "sld_id", all.x = T)
			sum_df$true_th21_around_ratio <- sum_df$true_around_th2/sum_df$true_around_th1
			sum_df$total_th21_ratio <- sum_df$total_th2/sum_df$total_th1
			sum_df$true_around_th2_pct <- sum_df$true_around_th2/sum_df$total_th2*100
			sum_df$true_around_th1_pct <- sum_df$true_around_th1/sum_df$total_th1*100
			sum_df$true_around_cd8t_pct <- sum_df$true_around_th2/sum_df$total_cd8t*100
			sum_df$true_around_th2_total_mdc <- sum_df$true_around_th2/nrow(fr_df)
			sum_df$true_around_th1_total_mdc <- sum_df$true_around_th1/nrow(fr_df)
			sum_df$true_around_cd8t_total_mdc <- sum_df$true_around_cd8t/nrow(fr_df)

			sum_df$total_th2_total_mdc <- sum_df$total_th2/nrow(fr_df)
			sum_df$total_th1_total_mdc <- sum_df$total_th1/nrow(fr_df)
			sum_df$total_cd8t_total_mdc <- sum_df$total_cd8t/nrow(fr_df)

			write.csv(sum_df, paste(fr_dir, "/merge_fr_", ifr, "_sample_wise_summary.csv", sep = ""))

			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_th21_around_ratio, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Th2/Th1 around mDC") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_true_th21_around_ratio.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = total_th21_ratio, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Th2/Th1 total") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_true_th21_total_ratio.png", sep = ""), dpi = 300, width = 3, height = 6)


			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_around_th2_pct, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Th2 percentage around mDC\n(To total Th2 cell number)") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_th2_pct_around_mdc.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_around_th1_pct, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Th1 percentage around mDC\n(To total Th1 cell number)") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_th1_pct_around_mdc.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_around_cd8t_pct, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "CD8T percentage around mDC\n(To total CD8T cell number)") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_cd8t_pct_around_mdc.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_around_th2_total_mdc, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Relative Th2 around mDC\n(To total mDC number)") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_th2_around_mdc_to_mdc.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_around_th1_total_mdc, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Relative Th1 around mDC\n(To total mDC number)") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_th1_around_mdc_to_mdc.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = true_around_cd8t_total_mdc, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Relative CD8T around mDC\n(To total mDC number)") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_cd8t_around_mdc_to_mdc.png", sep = ""), dpi = 300, width = 3, height = 6)

#			print(head(sum_df))
#			q(save = "no")
			box_gg <- ggplot(sum_df, aes(x = cohort, y = avg_around_th1, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Average Th1 cell numbers around mDC") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_avg_th1.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = avg_around_th2, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Average Th2 cell numbers around mDC") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_avg_th2.png", sep = ""), dpi = 300, width = 3, height = 6)

			box_gg <- ggplot(sum_df, aes(x = cohort, y = avg_around_cd8t, color = cohort)) +
				geom_boxplot() +
				geom_point(position = position_jitterdodge(dodge.width = 0.7)) +
				stat_compare_means(aes(group = cohort), label = "p.format") +
				labs(title = paste("Fixed radius:", ifr), y = "Average CD8+ T cell numbers around mDC") +
				theme_bw()
			ggsave(paste(fr_dir, "/merge_fr_", ifr, "_boxplot_avg_cd8t.png", sep = ""), dpi = 300, width = 3, height = 6)


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
			write.csv(gath_df, paste(fr_dir, "/merge_fr_", ifr, "_sample_wise_all_mdc_gather.csv", sep = ""))


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
			write.csv(gath_df, paste(fr_dir, "/merge_fr_", ifr, "_sample_wise_withthonly_mdc_gather.csv", sep = ""))
		}
	}
}

