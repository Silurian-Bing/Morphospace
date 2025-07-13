req_files <- c("OtoD_atrypidesB516.xlsx", "morphospace.xlsx")
if (!all(file.exists(req_files))) stop("Required diversity and morphospace data files missing")

rm(list = ls())
library(readxl);library(vegan);library(ggplot2);library(dplyr);library(forcats);library(mgcv);library(RColorBrewer)
library(viridis);library(tidyr);library(corrplot);library(gridExtra)
setwd("F:/HB/DataBase/Morphospace/")

calc_sqs <- function(occ_mtx, quota = 0.7, trials = 100) {
  if(nrow(occ_mtx) == 0 || ncol(occ_mtx) == 0) return(list(mean_div = 0, samples = numeric(0)))
  sp_counts <- colSums(occ_mtx); total_occ <- sum(sp_counts)
  if(total_occ == 0) return(list(mean_div = 0, samples = numeric(0)))
  sp_freq <- sp_counts / total_occ; target_spec <- floor(total_occ * quota)
  if(target_spec <= 0) return(list(mean_div = length(sp_counts[sp_counts > 0]), samples = rep(length(sp_counts[sp_counts > 0]), trials)))
  div_est <- numeric(trials)
  for(trial in 1:trials) {
    sampled_sp <- c(); spec_drawn <- 0
    while(spec_drawn < target_spec) {
      sp_idx <- sample(1:length(sp_freq), 1, prob = sp_freq)
      sampled_sp <- c(sampled_sp, sp_idx); spec_drawn <- spec_drawn + 1
    }
    div_est[trial] <- length(unique(sampled_sp))
  }
  return(list(mean_div = mean(div_est, na.rm = TRUE), samples = div_est[!is.na(div_est)]))
}

boot_chao1 <- function(dt_mtx, boots = 500) {
  n_sites <- nrow(dt_mtx)
  if (n_sites == 0 || ncol(dt_mtx) == 0) return(numeric(0))
  results <- numeric(boots)
  for(i in 1:boots) {
    boot_idx <- sample(1:n_sites, n_sites, replace = TRUE)
    boot_dt <- dt_mtx[boot_idx, , drop = FALSE]
    sp_occ <- colSums(boot_dt > 0); sp_occ <- sp_occ[sp_occ > 0]
    if(length(sp_occ) <= 1) { results[i] <- length(sp_occ); next }
    chao1_val <- tryCatch({ est_out <- vegan::estimateR(sp_occ); est_out[2] }, error = function(e) NA)
    results[i] <- chao1_val
  }
  return(results[!is.na(results)])
}

boot_morpho <- function(pca_scores, boots = 500) {
  if(nrow(pca_scores) < 3) return(list(volume = numeric(0), convex_hull = numeric(0), mean_dist = numeric(0)))
  volumes <- numeric(boots); convex_hulls <- numeric(boots); mean_dists <- numeric(boots)
  for(i in 1:boots) {
    boot_idx <- sample(1:nrow(pca_scores), nrow(pca_scores), replace = TRUE)
    boot_scores <- pca_scores[boot_idx, , drop = FALSE]
    if(ncol(boot_scores) >= 3) {
      pc_ranges <- apply(boot_scores[, 1:3], 2, function(x) max(x) - min(x))
      volumes[i] <- prod(pc_ranges)
    } else if(ncol(boot_scores) >= 2) {
      pc_ranges <- apply(boot_scores[, 1:2], 2, function(x) max(x) - min(x))
      volumes[i] <- prod(pc_ranges)
    } else { volumes[i] <- 0 }
    if(ncol(boot_scores) >= 2 && nrow(boot_scores) >= 3) {
      convex_hulls[i] <- tryCatch({
        hull_idx <- chull(boot_scores[, 1:2]); hull_pts <- boot_scores[hull_idx, 1:2]; n <- nrow(hull_pts)
        0.5 * abs(sum(hull_pts[1:n, 1] * hull_pts[c(2:n, 1), 2] - hull_pts[c(2:n, 1), 1] * hull_pts[1:n, 2]))
      }, error = function(e) 0)
    } else { convex_hulls[i] <- 0 }
    if(nrow(boot_scores) >= 2) {
      dist_mtx <- dist(boot_scores[, 1:min(2, ncol(boot_scores))]); mean_dists[i] <- mean(dist_mtx)
    } else { mean_dists[i] <- 0 }
  }
  return(list(volume = volumes[!is.na(volumes) & !is.infinite(volumes)],
              convex_hull = convex_hulls[!is.na(convex_hulls) & !is.infinite(convex_hulls)],
              mean_dist = mean_dists[!is.na(mean_dists) & !is.infinite(mean_dists)]))
}

div_file <- "OtoD_atrypidesB516.xlsx"; morph_file <- "morphospace.xlsx"
sht_names <- excel_sheets(div_file)
div_sum <- data.frame(Period = sht_names, TimeNumeric = 1:length(sht_names),
                      ObservedRichness = numeric(length(sht_names)), MeanChao1 = numeric(length(sht_names)),
                      MeanSQS = numeric(length(sht_names)), CombinedDiversity = numeric(length(sht_names)))

for (i in seq_along(sht_names)) {
  sht_name <- sht_names[i]
  curr_sht <- read_excel(div_file, sheet = sht_name, col_names = TRUE)
  if(nrow(curr_sht) == 0 || ncol(curr_sht) == 0) next
  
  df_calc <- as.data.frame(curr_sht)
  if (!is.numeric(df_calc[,1]) && !any(is.na(df_calc[,1])) && !any(duplicated(df_calc[,1]))) {
    rownames(df_calc) <- make.names(as.character(df_calc[,1]), unique = TRUE)
    df_calc <- df_calc[,-1, drop = FALSE]
  }
  
  dt_mtx_num <- suppressWarnings(apply(df_calc, 2, function(col) as.numeric(as.character(col))))
  if (is.vector(dt_mtx_num) && nrow(df_calc) == 1) {
    dt_mtx_num <- matrix(dt_mtx_num, nrow = 1, dimnames = list(rownames(df_calc)[1], colnames(df_calc)))
  } else if (!is.matrix(dt_mtx_num)){ dt_mtx_num <- as.matrix(dt_mtx_num); colnames(dt_mtx_num) <- colnames(df_calc) }
  dt_mtx_num[is.na(dt_mtx_num)] <- 0; dt_mtx_num[dt_mtx_num < 0] <- 0; dt_mtx_num[dt_mtx_num > 1] <- 1
  
  valid_loc <- rowSums(dt_mtx_num) > 0; valid_gen <- colSums(dt_mtx_num) > 0
  obs_rich <- sum(valid_gen); div_sum$ObservedRichness[i] <- obs_rich
  
  if(sum(valid_loc) > 0 && sum(valid_gen) > 1) {
    dt_analysis <- dt_mtx_num[valid_loc, valid_gen, drop = FALSE]
    chao1_est <- boot_chao1(dt_analysis, boots = 500)
    if(length(chao1_est) > 0) { div_sum$MeanChao1[i] <- mean(chao1_est, na.rm = TRUE)
    } else { div_sum$MeanChao1[i] <- obs_rich }
    
    sqs_res <- calc_sqs(dt_analysis, quota = 0.7, trials = 500)
    if(length(sqs_res$samples) > 0) { div_sum$MeanSQS[i] <- sqs_res$mean_div
    } else { div_sum$MeanSQS[i] <- obs_rich }
    
    if(length(chao1_est) > 0 && length(sqs_res$samples) > 0) {
      comb_samples <- c(chao1_est * 0.6, sqs_res$samples * 0.4)
      div_sum$CombinedDiversity[i] <- mean(comb_samples, na.rm = TRUE)
    } else { div_sum$CombinedDiversity[i] <- div_sum$MeanChao1[i] }
  } else { div_sum$MeanChao1[i] <- obs_rich; div_sum$MeanSQS[i] <- obs_rich; div_sum$CombinedDiversity[i] <- obs_rich }
}

morph_sht_names <- excel_sheets(morph_file)
morph_sum <- data.frame(Period = morph_sht_names, TimeNumeric = 1:length(morph_sht_names),
                        TotalVariance = numeric(length(morph_sht_names)), MorphoVolume = numeric(length(morph_sht_names)),
                        ConvexHullArea = numeric(length(morph_sht_names)), MeanPairwiseDistance = numeric(length(morph_sht_names)))

for (i in seq_along(morph_sht_names)) {
  sht_name <- morph_sht_names[i]
  tryCatch({
    df <- read_excel(path = morph_file, sheet = sht_name, range = cell_cols("A:AV"))
    if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) next
    
    df$TimePeriod <- sht_name; colnames(df)[1] <- "Genus"; df$Genus <- as.character(df$Genus)
    trait_cols <- names(df)[2:min(48, ncol(df)-1)]
    if (length(trait_cols) > 0) {
      df <- df %>% mutate(across(all_of(trait_cols), .fns = ~ {
        num_col <- suppressWarnings(as.numeric(as.character(.))); tidyr::replace_na(num_col, 0) }))
    }
    
    morph_dt <- df %>% select(all_of(trait_cols))
    if (any(is.na(morph_dt))) morph_dt[is.na(morph_dt)] <- 0
    zero_var_cols <- which(apply(morph_dt, 2, var) == 0)
    if (length(zero_var_cols) > 0) morph_dt <- morph_dt[, -zero_var_cols]
    if (nrow(morph_dt) < 3 || ncol(morph_dt) < 2) next
    
    pca_res <- rda(morph_dt, scale = TRUE); pca_sum <- summary(pca_res)
    exp_var <- pca_sum$cont$importance["Proportion Explained", ]
    morph_sum$TotalVariance[i] <- sum(exp_var[1:min(5, length(exp_var))]) * 100
    
    pca_scores <- scores(pca_res, display = "sites", choices = 1:min(3, ncol(morph_dt)))
    if(nrow(pca_scores) >= 3) {
      boot_res <- boot_morpho(pca_scores, boots = 500)
      if(length(boot_res$volume) > 0) {
        morph_sum$MorphoVolume[i] <- mean(boot_res$volume)
      } else {
        if(ncol(pca_scores) >= 3) { pc_ranges <- apply(pca_scores[, 1:3], 2, function(x) max(x) - min(x))
        } else if(ncol(pca_scores) >= 2) { pc_ranges <- apply(pca_scores[, 1:2], 2, function(x) max(x) - min(x)) }
        morph_sum$MorphoVolume[i] <- prod(pc_ranges)
      }
      if(length(boot_res$mean_dist) > 0) { morph_sum$MeanPairwiseDistance[i] <- mean(boot_res$mean_dist)
      } else { dist_mtx <- dist(pca_scores[, 1:min(2, ncol(pca_scores))]); morph_sum$MeanPairwiseDistance[i] <- mean(dist_mtx) }
      if(length(boot_res$convex_hull) > 0) { morph_sum$ConvexHullArea[i] <- mean(boot_res$convex_hull)
      } else { morph_sum$ConvexHullArea[i] <- tryCatch({
        hull_idx <- chull(pca_scores[, 1:2]); hull_pts <- pca_scores[hull_idx, 1:2]; n <- nrow(hull_pts)
        0.5 * abs(sum(hull_pts[1:n, 1] * hull_pts[c(2:n, 1), 2] - hull_pts[c(2:n, 1), 1] * hull_pts[1:n, 2]))
      }, error = function(e) 0) }
    }
  }, error = function(e) {})
}

common_periods <- intersect(div_sum$Period, morph_sum$Period)
if(length(common_periods) == 0) stop("No common time periods found")

integ_dt <- merge(div_sum[div_sum$Period %in% common_periods, ], morph_sum[morph_sum$Period %in% common_periods, ],
                  by = "Period", suffixes = c("_div", "_morpho")) %>% arrange(TimeNumeric_div) %>%
  mutate(Period = factor(Period, levels = Period, ordered = TRUE), TimeNumeric = row_number())

gam_dt <- integ_dt %>% mutate(Div_Chao1_scaled = scale(MeanChao1)[,1], Div_SQS_scaled = scale(MeanSQS)[,1],
                              Div_Comb_scaled = scale(CombinedDiversity)[,1], MorphoVol_scaled = scale(MorphoVolume)[,1])

gam_div_chao1 <- gam(MeanChao1 ~ s(TimeNumeric, k = min(4, nrow(gam_dt)-1)), data = gam_dt, method = "REML")
gam_div_sqs <- gam(MeanSQS ~ s(TimeNumeric, k = min(4, nrow(gam_dt)-1)), data = gam_dt, method = "REML")
gam_div_comb <- gam(CombinedDiversity ~ s(TimeNumeric, k = min(4, nrow(gam_dt)-1)), data = gam_dt, method = "REML")
gam_morpho_vol <- gam(MorphoVolume ~ s(TimeNumeric, k = min(4, nrow(gam_dt)-1)), data = gam_dt, method = "REML")

create_div_comp_plot <- function() {
  chao1_col <- "#1f77b4"; sqs_col <- "#ff7f0e"; comb_col <- "#2ca02c"
    
  ggplot(integ_dt, aes(x = TimeNumeric)) +
    geom_line(aes(y = MeanChao1), color = chao1_col, linewidth = 1.5, alpha = 0.8) +
    geom_point(aes(y = MeanChao1), color = chao1_col, size = 3.5, alpha = 0.8) +
    geom_line(aes(y = MeanSQS), color = sqs_col, linewidth = 1.5, alpha = 0.8) +
    geom_point(aes(y = MeanSQS), color = sqs_col, size = 3.5, shape = 17, alpha = 0.8) +
    geom_line(aes(y = CombinedDiversity), color = comb_col, linewidth = 1.5, alpha = 0.8) +
    geom_point(aes(y = CombinedDiversity), color = comb_col, size = 3.5, shape = 15, alpha = 0.8) +
    scale_x_continuous(breaks = 1:nrow(integ_dt), labels = integ_dt$Period) +
    labs(title = "Diversity Estimates Comparison: Chao1, SQS, and Combined Methods",
         subtitle = "Circle: Chao1 • Triangle: SQS • Square: Combined", x = "Geological Time Period", y = "Diversity Estimate") +
    theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                                          axis.text.x = element_text(angle = 45, hjust = 1))
}

create_dual_axis_plot <- function() {
  p_dt <- integ_dt %>% mutate(scale_factor = max(CombinedDiversity, na.rm = TRUE) / max(MorphoVolume, na.rm = TRUE),
                              MorphoVol_scaled = MorphoVolume * scale_factor)
  div_col <- "#2ca02c"; morpho_col <- "#ff7f0e"
    
  ggplot(p_dt, aes(x = TimeNumeric)) +
    geom_line(aes(y = CombinedDiversity), color = div_col, linewidth = 1.8, alpha = 0.9) +
    geom_point(aes(y = CombinedDiversity), color = div_col, size = 4.5, alpha = 0.8) +
    geom_line(aes(y = MorphoVol_scaled), color = morpho_col, linewidth = 1.8, alpha = 0.9) +
    geom_point(aes(y = MorphoVol_scaled), color = morpho_col, size = 4.5, shape = 17, alpha = 0.8) +
    scale_x_continuous(breaks = 1:nrow(integ_dt), labels = integ_dt$Period) +
    scale_y_continuous(name = "Combined Diversity", sec.axis = sec_axis(~ . / p_dt$scale_factor[1], name = "Morphospace Volume")) +
    labs(title = "Combined Diversity vs Morphospace Occupancy Through Time", x = "Geological Time Period") +
    theme_minimal(base_size = 14) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                                          axis.title.y.left = element_text(color = div_col, face = "bold"),
                                          axis.title.y.right = element_text(color = morpho_col, face = "bold"),
                                          axis.text.x = element_text(angle = 45, hjust = 1))
}

p1 <- create_div_comp_plot(); p2 <- create_dual_axis_plot()
ggsave("diversity_methods_comparison.png", p1, width = 16, height = 10, dpi = 300, bg = "white")
ggsave("combined_diversity_vs_morphospace.png", p2, width = 16, height = 10, dpi = 300, bg = "white")

corr_res <- cor(integ_dt[, c("MeanChao1", "MeanSQS", "CombinedDiversity", "MorphoVolume", "ConvexHullArea", "MeanPairwiseDistance")], 
                use = "complete.obs")

write.csv(integ_dt, "integrated_diversity_morphospace_data.csv", row.names = FALSE)

sink("morphospace_diversity_analysis_results.txt")
cat("Morphospace vs Diversity Analysis Results\n===========================================\n\n")
cat("Common time periods:", length(common_periods), "\nPeriods:", paste(common_periods, collapse = ", "), "\n\n")
cat("Diversity Statistics:\nChao1:\n"); print(summary(integ_dt$MeanChao1))
cat("\nSQS:\n"); print(summary(integ_dt$MeanSQS))
cat("\nCombined:\n"); print(summary(integ_dt$CombinedDiversity))
cat("\nMorphospace Statistics:\nVolume:\n"); print(summary(integ_dt$MorphoVolume))
cat("\nCorrelation Matrix:\n"); print(round(corr_res, 3))
sink()