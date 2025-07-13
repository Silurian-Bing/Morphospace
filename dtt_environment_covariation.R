req_files <- c("tem470_370.csv", "sea470_370.csv", "dtt_timebin_data.csv")
if (!all(file.exists(req_files))) stop("Required environmental data files missing")

library(readr);library(dplyr);library(tidyr);library(ggplot2);library(viridis);library(mgcv);library(gridExtra);library(scales);library(RColorBrewer)
setwd("F:/HB/DataBase/Morphospace/"); options(scipen = 999)

gam_smoother <- function(dt, age_col, val_col, grid, k_val = 15) {
  dt_flt <- dt[dt[[age_col]] >= min(grid) & dt[[age_col]] <= max(grid), ]
  dt_flt <- dt_flt[!is.na(dt_flt[[val_col]]), ]
  if (nrow(dt_flt) < 4) return(list(fit = approx(dt_flt[[age_col]], dt_flt[[val_col]], xout = grid)$y, se = NULL))
  gam_dt <- data.frame(age = dt_flt[[age_col]], value = dt_flt[[val_col]])
  k_use <- min(k_val, max(4, nrow(gam_dt) - 2))
  gam_mod <- gam(value ~ s(age, k = k_use, bs = "cs"), data = gam_dt, method = "REML")
  pred_se <- predict(gam_mod, newdata = data.frame(age = grid), se.fit = TRUE)
  return(list(fit = pred_se$fit, se = pred_se$se.fit, lower = pred_se$fit - 1.96 * pred_se$se.fit, upper = pred_se$fit + 1.96 * pred_se$se.fit))
}

norm_01 <- function(x) { if (all(is.na(x)) || length(unique(x[!is.na(x)])) <= 1) return(x); (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) }

temp_dt <- read_csv("tem470_370.csv"); colnames(temp_dt) <- c("Age", "GAT", "Tropical")
sea_dt <- read_csv("sea470_370.csv", locale = locale(encoding = "UTF-8")); colnames(sea_dt) <- c("Age", "Sealevel")

if (file.exists("dtt_timebin_data.csv")) {
  dtt_dt <- read_csv("dtt_timebin_data.csv")
} else { stop("DTT data file not found") }

comm_min_age <- max(min(temp_dt$Age), min(sea_dt$Age), min(dtt_dt$time_midpoint))
comm_max_age <- min(max(temp_dt$Age), max(sea_dt$Age), max(dtt_dt$time_midpoint))
time_grid <- seq(comm_min_age, comm_max_age, by = 0.5)

gat_smooth <- gam_smoother(temp_dt, "Age", "GAT", time_grid)
trop_smooth <- gam_smoother(temp_dt, "Age", "Tropical", time_grid)
sea_smooth <- gam_smoother(sea_dt, "Age", "Sealevel", time_grid)

dtt_flt <- dtt_dt[dtt_dt$time_midpoint >= comm_min_age & dtt_dt$time_midpoint <= comm_max_age, ]
excl_vars <- c("time_midpoint", "time_start", "time_end", "time_bin", "pc1_range", "pc2_range", "mean_pairwise")
all_vars <- names(dtt_flt); dtt_vars <- all_vars[!all_vars %in% excl_vars]
num_vars <- sapply(dtt_flt[dtt_vars], is.numeric); dtt_vars <- dtt_vars[num_vars]

dtt_smoothed <- list()
for (var in dtt_vars) {
  var_dt <- dtt_flt[[var]]
  if (all(is.na(var_dt)) || length(unique(var_dt[!is.na(var_dt)])) <= 1) next
  smooth_res <- gam_smoother(dtt_flt, "time_midpoint", var, time_grid)
  dtt_smoothed[[var]] <- smooth_res$fit
}

complete_dt <- data.frame(age = time_grid, Global_Avg_Temp = gat_smooth$fit, Tropical_Temp = trop_smooth$fit, Sea_Level = sea_smooth$fit)
for (var in names(dtt_smoothed)) { complete_dt[[var]] <- dtt_smoothed[[var]] }
complete_dt_clean <- complete_dt[complete.cases(complete_dt), ]

var_mapping <- list("Global_Avg_Temp" = "Global Avg Temp", "Tropical_Temp" = "Tropical Temp", "Sea_Level" = "Sea Level",
                    "total_var" = "Total Disparity", "pc1_var" = "PC1 Disparity", "pc2_var" = "PC2 Disparity", 
                    "n_genera" = "Diversity", "diversity" = "Diversity")

corr_data <- complete_dt_clean[, names(complete_dt_clean) != "age"]
final_names <- names(corr_data)
for (i in seq_along(final_names)) { if (final_names[i] %in% names(var_mapping)) final_names[i] <- var_mapping[[final_names[i]]] }
names(corr_data) <- final_names

full_cor_mtx <- cor(corr_data, use = "complete.obs", method = "pearson")

morph_var <- NULL; div_var <- NULL
morph_cands <- c("total_var", "total_disparity", "morphospace", "variance")
for (cand in morph_cands) { if (cand %in% names(dtt_smoothed)) { morph_var <- cand; break } }
div_cands <- c("n_genera", "diversity", "richness", "n_taxa")
for (cand in div_cands) { if (cand %in% names(dtt_smoothed)) { div_var <- cand; break } }

if (is.null(morph_var) && length(dtt_smoothed) > 0) morph_var <- names(dtt_smoothed)[1]
if (is.null(div_var) && length(dtt_smoothed) > 1) div_var <- names(dtt_smoothed)[2]
else if (is.null(div_var)) div_var <- morph_var

gat_res <- gam_smoother(temp_dt, "Age", "GAT", time_grid)
sea_res <- gam_smoother(sea_dt, "Age", "Sealevel", time_grid)
morph_res <- if(!is.null(morph_var)) gam_smoother(dtt_flt, "time_midpoint", morph_var, time_grid) else list(fit = rep(0, length(time_grid)))
div_res <- if(!is.null(div_var)) gam_smoother(dtt_flt, "time_midpoint", div_var, time_grid) else list(fit = rep(0, length(time_grid)))

core_dt <- data.frame(age = time_grid, temperature = gat_res$fit, sealevel = sea_res$fit, 
                      morphospace = morph_res$fit, diversity = div_res$fit,
                      temp_lower = if(!is.null(gat_res$lower)) gat_res$lower else gat_res$fit,
                      temp_upper = if(!is.null(gat_res$upper)) gat_res$upper else gat_res$fit,
                      sea_lower = if(!is.null(sea_res$lower)) sea_res$lower else sea_res$fit,
                      sea_upper = if(!is.null(sea_res$upper)) sea_res$upper else sea_res$fit)

core_dt_comp <- core_dt[complete.cases(core_dt[, c("age", "temperature", "sealevel", "morphospace", "diversity")]), ]

core_vars <- c("temperature", "sealevel", "morphospace", "diversity")
cor_mtx <- cor(core_dt_comp[core_vars], use = "complete.obs", method = "pearson")

vibrant_cols <- c("Global Average Temperature" = "#E31A1C", "Sea Level" = "#1F78B4", 
                  "Taxonomic Diversity" = "#33A02C", "Morphospace Volume" = "#FF7F00")

core_dt_norm <- core_dt_comp; core_dt_norm[, core_vars] <- lapply(core_dt_comp[, core_vars], norm_01)

var_labels <- c("temperature" = "Global Average Temperature", "sealevel" = "Sea Level",
                "diversity" = "Taxonomic Diversity", "morphospace" = "Morphospace Volume")

plot_dt <- data.frame()
for (var in core_vars) {
  var_dt <- data.frame(age = core_dt_norm$age, variable = var, variable_label = var_labels[var],
                       value = core_dt_norm[[var]], lower = core_dt_norm[[var]], upper = core_dt_norm[[var]])
  plot_dt <- rbind(plot_dt, var_dt)
}

plot_dt$variable_label <- factor(plot_dt$variable_label, levels = names(vibrant_cols))

theme_prof <- theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 14, face = "bold"),
        panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", fill = NA))

main_plot <- ggplot(plot_dt, aes(x = age, y = value, color = variable_label, fill = variable_label)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
  geom_line(size = 2.5, alpha = 0.9) + geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = vibrant_cols, name = "Variables") + scale_fill_manual(values = vibrant_cols, guide = "none") +
  scale_x_reverse(name = "Age (Ma)", breaks = pretty_breaks(n = 8)) +
  scale_y_continuous(name = "Normalized Values", breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_prof + labs(title = "DTT-Environment Covariation Analysis", 
                    subtitle = "Four Key Variables Through Time with GAM Smoothing") +
  guides(color = guide_legend(title = "Variables", override.aes = list(size = 4)))

create_cor_heatmap <- function(cor_mtx) {
  cor_dt <- expand.grid(Var1 = rownames(cor_mtx), Var2 = colnames(cor_mtx))
  cor_dt$value <- as.vector(cor_mtx)
  cor_dt$Var1 <- factor(cor_dt$Var1, levels = rownames(cor_mtx))
  cor_dt$Var2 <- factor(cor_dt$Var2, levels = rev(rownames(cor_mtx)))
  
  ggplot(cor_dt, aes(x = Var1, y = Var2, fill = value)) + geom_tile(color = "white", size = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value)), color = "white", size = 4.5, fontface = "bold") +
    scale_fill_gradient2(low = "#3F4F96", mid = "white", high = "#C73E1D", midpoint = 0, limit = c(-1, 1),
                         name = "Correlation\nCoefficient") +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
                                          panel.grid = element_blank(), axis.title = element_blank()) +
    labs(title = "Variable Correlation Matrix") + coord_fixed(ratio = 1)
}

cor_plot <- create_cor_heatmap(full_cor_mtx)

ggsave("dtt_environment_comprehensive.png", main_plot, width = 14, height = 8, dpi = 300, bg = "white")
ggsave("dtt_environment_correlation.png", cor_plot, width = 10, height = 8, dpi = 300, bg = "white")

write_csv(complete_dt_clean, "dtt_environment_complete_data.csv")
write_csv(core_dt_comp, "dtt_environment_core_data.csv")
write.csv(full_cor_mtx, "dtt_environment_correlation_matrix.csv")