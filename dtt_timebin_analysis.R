req_files <- c("morphospace.xlsx", "taxon_age_data.csv")
if (!all(file.exists(req_files))) stop("Required morphospace data files missing")

library(readxl);library(dplyr);library(tidyr);library(ggplot2);library(viridis);library(vegan);library(gridExtra)
setwd("F:/HB/DataBase/Morphospace/")

assign_time_bins <- function(fad, lad, bins) {
  present <- c()
  for (i in 1:(length(bins)-1)) {
    bin_start <- bins[i+1]; bin_end <- bins[i]
    if (fad >= bin_end && lad <= bin_start) present <- c(present, paste0(bin_start, "-", bin_end, "Ma"))
  }
  return(present)
}

calc_disparity_metrics <- function(genera, morph_dt) {
  if (length(genera) < 2) return(list(n_genera = length(genera), pc1_var = NA, pc2_var = NA, total_var = NA, mean_pairwise = NA))
  avail <- intersect(genera, morph_dt$Genus)
  if (length(avail) < 2) return(list(n_genera = length(avail), pc1_var = NA, pc2_var = NA, total_var = NA, mean_pairwise = NA))
  subset_dt <- morph_dt[morph_dt$Genus %in% avail, ]
  pc1_var <- var(subset_dt$PC1_mean, na.rm = TRUE); pc2_var <- var(subset_dt$PC2_mean, na.rm = TRUE)
  pc_dt <- subset_dt[, c("PC1_mean", "PC2_mean")]; pc_dt <- pc_dt[complete.cases(pc_dt), ]
  if (nrow(pc_dt) >= 2) {
    total_var <- sum(apply(pc_dt, 2, var, na.rm = TRUE), na.rm = TRUE)
    mean_pairwise <- mean(dist(pc_dt), na.rm = TRUE)
  } else { total_var <- NA; mean_pairwise <- NA }
  return(list(n_genera = length(avail), pc1_var = pc1_var, pc2_var = pc2_var, total_var = total_var, mean_pairwise = mean_pairwise))
}

file_path <- "morphospace.xlsx"; sht_names <- excel_sheets(file_path)
all_dt_list <- lapply(sht_names, function(sht) {
  df <- tryCatch({ read_excel(path = file_path, sheet = sht, range = cell_cols("A:AV")) }, 
                 error = function(e) return(NULL))
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$TimePeriod <- sht; colnames(df)[1] <- "Genus"; df$Genus <- as.character(df$Genus)
  trait_cols <- names(df)[2:48]
  if (length(trait_cols) > 0) {
    df <- df %>% mutate(across(all_of(trait_cols), .fns = ~ {
      num_col <- suppressWarnings(as.numeric(as.character(.)))
      tidyr::replace_na(num_col, 0) }))
  }
  return(df)
})

all_dt_list <- all_dt_list[!sapply(all_dt_list, is.null)]; comb_dt <- bind_rows(all_dt_list)
trait_cols <- names(all_dt_list[[1]])[2:48]; char_cols <- intersect(trait_cols, names(comb_dt))
morph_dt <- comb_dt %>% select(all_of(char_cols)); meta_dt <- comb_dt %>% select(Genus, TimePeriod)

if (any(is.na(morph_dt))) morph_dt[is.na(morph_dt)] <- 0
zero_var_cols <- which(apply(morph_dt, 2, var) == 0)
if (length(zero_var_cols) > 0) morph_dt <- morph_dt[, -zero_var_cols]

pca_res <- prcomp(morph_dt, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_res$x[, 1:min(5, ncol(pca_res$x))])
exp_var <- summary(pca_res)$importance["Proportion of Variance", ]

plot_dt <- bind_cols(meta_dt, pca_scores)
plot_dt$TimePeriod <- factor(plot_dt$TimePeriod, levels = sht_names, ordered = TRUE)

rng <- read.csv("taxon_age_data.csv"); min_age <- min(rng$LAD, na.rm = TRUE); max_age <- max(rng$FAD, na.rm = TRUE)
bin_size <- 5; time_bins <- seq(from = floor(min_age/bin_size) * bin_size, to = ceiling(max_age/bin_size) * bin_size, by = bin_size)

genus_time_assign <- list()
for (i in 1:nrow(rng)) {
  genus <- rng$taxon[i]; fad <- rng$FAD[i]; lad <- rng$LAD[i]
  if (!is.na(fad) && !is.na(lad)) {
    assigned_bins <- assign_time_bins(fad, lad, time_bins)
    if (length(assigned_bins) > 0) {
      for (bin in assigned_bins) {
        if (is.null(genus_time_assign[[bin]])) genus_time_assign[[bin]] <- c()
        genus_time_assign[[bin]] <- c(genus_time_assign[[bin]], genus)
      }
    }
  }
}

genus_morph_sum <- plot_dt %>% group_by(Genus) %>%
  summarise(PC1_mean = mean(PC1, na.rm = TRUE), PC2_mean = mean(PC2, na.rm = TRUE),
            PC3_mean = if("PC3" %in% names(.)) mean(PC3, na.rm = TRUE) else NA, .groups = 'drop')

disp_res <- list(); valid_bins <- names(genus_time_assign)[sapply(genus_time_assign, length) >= 2]
for (bin_name in valid_bins) {
  genera_in_bin <- genus_time_assign[[bin_name]]
  disparity <- calc_disparity_metrics(genera_in_bin, genus_morph_sum)
  bin_ages <- as.numeric(strsplit(gsub("Ma", "", bin_name), "-")[[1]])
  bin_midpoint <- mean(bin_ages)
  disp_res[[bin_name]] <- c(time_midpoint = bin_midpoint, time_bin = bin_name, disparity)
}

dtt_dt <- do.call(rbind, lapply(disp_res, function(x) {
  data.frame(time_midpoint = as.numeric(x$time_midpoint), time_bin = as.character(x$time_bin),
             n_genera = as.numeric(x$n_genera), pc1_var = as.numeric(x$pc1_var), pc2_var = as.numeric(x$pc2_var),
             total_var = as.numeric(x$total_var), mean_pairwise = as.numeric(x$mean_pairwise), stringsAsFactors = FALSE)
}))

dtt_dt <- dtt_dt[order(dtt_dt$time_midpoint, decreasing = TRUE), ]

create_dtt_plots <- function(dt) {
  p1 <- ggplot(dt, aes(x = time_midpoint, y = total_var)) + geom_line(color = "blue", size = 1.2) +
    geom_point(color = "darkblue", size = 3) + scale_x_reverse() +
    labs(title = "Disparity Through Time - Total Variance", x = "Time (Ma)", y = "Total Variance") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p2 <- ggplot(dt, aes(x = time_midpoint, y = pc1_var)) + geom_line(color = "red", size = 1.2) +
    geom_point(color = "darkred", size = 3) + scale_x_reverse() +
    labs(title = "DTT - PC1 Variance", x = "Time (Ma)", y = "PC1 Variance") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p3 <- ggplot(dt, aes(x = time_midpoint, y = pc2_var)) + geom_line(color = "green", size = 1.2) +
    geom_point(color = "darkgreen", size = 3) + scale_x_reverse() +
    labs(title = "DTT - PC2 Variance", x = "Time (Ma)", y = "PC2 Variance") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p4 <- ggplot(dt, aes(x = time_midpoint, y = n_genera)) + geom_line(color = "orange", size = 1.2) +
    geom_point(color = "darkorange", size = 3) + scale_x_reverse() +
    labs(title = "Taxonomic Diversity Through Time", x = "Time (Ma)", y = "Number of Genera") +
    theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(list(p1, p2, p3, p4))
}

plots <- create_dtt_plots(dtt_dt)
ggsave("dtt_total_variance.png", plots[[1]], width = 10, height = 6, dpi = 300)
ggsave("dtt_pc1_variance.png", plots[[2]], width = 10, height = 6, dpi = 300)
ggsave("dtt_pc2_variance.png", plots[[3]], width = 10, height = 6, dpi = 300)
ggsave("diversity_through_time.png", plots[[4]], width = 10, height = 6, dpi = 300)

png("dtt_comprehensive_timebins.png", width = 1600, height = 1200, res = 300)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2)
dev.off()

analyze_trends <- function(dt) {
  if (nrow(dt) < 3) return(data.frame(Metric = "Insufficient_data", Trend = "Cannot_analyze", P_value = NA))
  
  total_lm <- lm(total_var ~ time_midpoint, data = dt, na.action = na.omit)
  pc1_lm <- lm(pc1_var ~ time_midpoint, data = dt, na.action = na.omit)
  div_lm <- lm(n_genera ~ time_midpoint, data = dt)
  
  get_trend <- function(model) {
    p_val <- summary(model)$coefficients[2, 4]
    if (p_val < 0.05) return(ifelse(coef(model)[2] > 0, "Increasing", "Decreasing"))
    return("No significant trend")
  }
  
  return(data.frame(Metric = c("Total_Variance", "PC1_Variance", "Diversity"),
                    Trend = c(get_trend(total_lm), get_trend(pc1_lm), get_trend(div_lm)),
                    P_value = c(summary(total_lm)$coefficients[2, 4], summary(pc1_lm)$coefficients[2, 4], summary(div_lm)$coefficients[2, 4])))
}

trend_analysis <- analyze_trends(dtt_dt)

write.csv(dtt_dt, "dtt_timebin_data.csv", row.names = FALSE)
write.csv(trend_analysis, "dtt_trend_analysis.csv", row.names = FALSE)

dtt_summary <- list(dtt_data = dtt_dt, trend_analysis = trend_analysis, time_bins = time_bins,
                    bin_size = bin_size, total_time_range = c(min_age, max_age))
saveRDS(dtt_summary, "dtt_timebin_analysis.rds")