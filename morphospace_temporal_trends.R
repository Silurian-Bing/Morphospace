library(readxl);library(dplyr);library(ggplot2);library(gridExtra);library(mgcv);library(viridis);library(corrplot);library(changepoint)
setwd("F:/HB/DataBase/Morphospace/")

age_data <- read.csv("taxon_age_data.csv")
excel_file <- "morphospace.xlsx"; sheet_names <- excel_sheets(excel_file)
all_morpho_data <- list(); all_taxa <- c()

for (sheet in sheet_names) {
  current_data <- read_excel(excel_file, sheet = sheet)
  if (ncol(current_data) > 1) {
    taxa_col <- current_data[[1]]; valid_rows <- !is.na(taxa_col) & taxa_col != ""
    if (sum(valid_rows) > 0) {
      sheet_data <- current_data[valid_rows, ]; names(sheet_data)[1] <- "taxon"
      sheet_data$geological_stage <- sheet; all_morpho_data[[sheet]] <- sheet_data
      all_taxa <- c(all_taxa, sheet_data$taxon)
    }
  }
}

common_taxa <- intersect(age_data$taxon, unique(all_taxa))
if (length(common_taxa) < 10) stop("Insufficient taxa for analysis!")

if (length(all_morpho_data) > 0) {
  all_colnames <- lapply(all_morpho_data, names); common_columns <- Reduce(intersect, all_colnames)
  final_morpho_data <- data.frame()
  for (sheet_data in all_morpho_data) {
    sheet_subset <- sheet_data[, common_columns]; final_morpho_data <- rbind(final_morpho_data, sheet_subset)
  }
  final_morpho_data <- final_morpho_data[final_morpho_data$taxon %in% common_taxa, ]
}

morpho_features <- final_morpho_data[, !names(final_morpho_data) %in% c("taxon", "geological_stage")]
morpho_matrix <- matrix(NA, nrow = nrow(morpho_features), ncol = ncol(morpho_features))
colnames(morpho_matrix) <- names(morpho_features); rownames(morpho_matrix) <- final_morpho_data$taxon

valid_cols <- logical(ncol(morpho_features))
for (i in 1:ncol(morpho_features)) {
  col_data <- suppressWarnings(as.numeric(morpho_features[[i]]))
  if (sum(!is.na(col_data)) >= 3 && var(col_data, na.rm = TRUE) > 0) {
    morpho_matrix[, i] <- col_data; valid_cols[i] <- TRUE
  }
}

morpho_matrix <- morpho_matrix[, valid_cols, drop = FALSE]
valid_rows <- rowSums(!is.na(morpho_matrix)) >= 3
morpho_matrix <- morpho_matrix[valid_rows, , drop = FALSE]
available_taxa <- rownames(morpho_matrix); final_taxa <- intersect(available_taxa, age_data$taxon)
final_morpho_clean <- morpho_matrix[final_taxa, , drop = FALSE]

for (i in 1:ncol(final_morpho_clean)) {
  na_indices <- is.na(final_morpho_clean[, i])
  if (any(na_indices)) {
    col_mean <- mean(final_morpho_clean[, i], na.rm = TRUE)
    final_morpho_clean[na_indices, i] <- col_mean
  }
}

final_age_data <- age_data[age_data$taxon %in% final_taxa, ]
final_age_data <- final_age_data[match(final_taxa, final_age_data$taxon), ]
final_age_data$avg_age <- (final_age_data$FAD + final_age_data$LAD) / 2

# PCO (Principal Coordinates Analysis)
dist_matrix <- dist(scale(final_morpho_clean, center = TRUE, scale = TRUE))
pco_result <- cmdscale(dist_matrix, k = min(nrow(final_morpho_clean) - 1, ncol(final_morpho_clean)), eig = TRUE)
variance_explained <- (pco_result$eig / sum(pco_result$eig)) * 100
n_pcos <- min(6, ncol(pco_result$points)); pco_scores <- pco_result$points[, 1:n_pcos]

morpho_data <- data.frame(
  taxon = rownames(pco_scores), PCO1 = pco_scores[, 1],
  PCO2 = if(n_pcos > 1) pco_scores[, 2] else 0, PCO3 = if(n_pcos > 2) pco_scores[, 3] else 0,
  PCO4 = if(n_pcos > 3) pco_scores[, 4] else 0, PCO5 = if(n_pcos > 4) pco_scores[, 5] else 0,
  PCO6 = if(n_pcos > 5) pco_scores[, 6] else 0, avg_age = final_age_data$avg_age,
  FAD = final_age_data$FAD, LAD = final_age_data$LAD,
  age_range = final_age_data$LAD - final_age_data$FAD, stringsAsFactors = FALSE
)

morpho_data$morpho_space_volume_2d <- abs(morpho_data$PCO1 * morpho_data$PCO2)
morpho_data$morpho_space_volume_3d <- abs(morpho_data$PCO1 * morpho_data$PCO2 * morpho_data$PCO3)
pco_center <- apply(morpho_data[, c("PCO1", "PCO2", "PCO3")], 2, mean)
morpho_data$distance_to_center <- sqrt(rowSums((morpho_data[, c("PCO1", "PCO2", "PCO3")] - 
                                                  matrix(pco_center, nrow = nrow(morpho_data), ncol = 3, byrow = TRUE))^2))
morpho_data$morpho_disparity <- sqrt(rowSums(morpho_data[, paste0("PCO", 1:min(4, n_pcos))]^2))

morpho_data <- morpho_data[order(morpho_data$avg_age, decreasing = TRUE), ]
morpho_data$pco1_rate <- c(NA, diff(morpho_data$PCO1) / diff(morpho_data$avg_age))
morpho_data$pco2_rate <- c(NA, diff(morpho_data$PCO2) / diff(morpho_data$avg_age))
morpho_data$overall_rate <- sqrt(morpho_data$pco1_rate^2 + morpho_data$pco2_rate^2)

pco_weights <- variance_explained[1:min(4, n_pcos)] / sum(variance_explained[1:min(4, n_pcos)])
morpho_data$weighted_complexity <- 0
for (i in 1:length(pco_weights)) {
  morpho_data$weighted_complexity <- morpho_data$weighted_complexity + 
    abs(morpho_data[[paste0("PCO", i)]]) * pco_weights[i]
}
morpho_data$specialization <- apply(abs(morpho_data[, paste0("PCO", 1:min(3, n_pcos))]), 1, max)

time_windows <- seq(from = ceiling(max(morpho_data$avg_age)), to = floor(min(morpho_data$avg_age)), by = -10)
morpho_data$time_window <- cut(morpho_data$avg_age, breaks = rev(time_windows), 
                               labels = paste0(time_windows[-length(time_windows)], "-", time_windows[-1], "Ma"),
                               include.lowest = TRUE)

trend_vars <- c("PCO1", "PCO2", "PCO3", "morpho_disparity", "distance_to_center", "weighted_complexity", "specialization")
trend_results <- data.frame(variable = trend_vars, linear_r = NA, linear_p = NA, kendall_tau = NA, 
                            kendall_p = NA, direction = NA, strength = NA, stringsAsFactors = FALSE)

for (i in 1:length(trend_vars)) {
  var_name <- trend_vars[i]; var_data <- morpho_data[[var_name]]
  linear_test <- cor.test(morpho_data$avg_age, var_data, method = "pearson")
  trend_results$linear_r[i] <- linear_test$estimate; trend_results$linear_p[i] <- linear_test$p.value
  kendall_test <- cor.test(morpho_data$avg_age, var_data, method = "kendall")
  trend_results$kendall_tau[i] <- kendall_test$estimate; trend_results$kendall_p[i] <- kendall_test$p.value
  trend_results$direction[i] <- ifelse(linear_test$estimate > 0, "positive", "negative")
  trend_results$strength[i] <- ifelse(abs(linear_test$estimate) > 0.5, "strong",
                                      ifelse(abs(linear_test$estimate) > 0.3, "moderate", "weak"))
}

gam_results <- list()
gam_summary <- data.frame(variable = trend_vars, gam_r2 = NA, gam_p = NA, linear_r2 = NA, 
                          nonlinear_improvement = NA, stringsAsFactors = FALSE)

for (i in 1:length(trend_vars)) {
  var_name <- trend_vars[i]; var_data <- morpho_data[[var_name]]
  gam_model <- gam(var_data ~ s(avg_age, k = 5), data = morpho_data); gam_results[[var_name]] <- gam_model
  linear_model <- lm(var_data ~ avg_age, data = morpho_data)
  gam_summary$gam_r2[i] <- summary(gam_model)$r.sq; gam_summary$gam_p[i] <- summary(gam_model)$s.pv
  gam_summary$linear_r2[i] <- summary(linear_model)$r.squared
  gam_summary$nonlinear_improvement[i] <- gam_summary$gam_r2[i] - gam_summary$linear_r2[i]
  morpho_data[[paste0(var_name, "_gam_fit")]] <- predict(gam_model)
}

changepoint_results <- list()
for (var_name in c("PCO1", "PCO2", "morpho_disparity")) {
  var_data <- morpho_data[[var_name]]
  tryCatch({
    cpt_mean <- cpt.mean(var_data, method = "PELT"); changepoints <- cpts(cpt_mean)
    if (length(changepoints) > 0) {
      change_ages <- morpho_data$avg_age[changepoints]; changepoint_results[[var_name]] <- change_ages
    }
  }, error = function(e) {})
}

window_stats <- morpho_data %>% filter(!is.na(time_window)) %>% group_by(time_window) %>%
  summarise(n_taxa = n(), avg_age_mean = mean(avg_age), pco1_mean = mean(PCO1), pco1_sd = sd(PCO1),
            pco2_mean = mean(PCO2), pco2_sd = sd(PCO2), morpho_disparity_mean = mean(morpho_disparity),
            morpho_disparity_sd = sd(morpho_disparity), space_volume_mean = mean(morpho_space_volume_2d),
            specialization_mean = mean(specialization), .groups = 'drop') %>% arrange(desc(avg_age_mean))

window_stats$pco1_change_rate <- c(NA, diff(window_stats$pco1_mean) / diff(window_stats$avg_age_mean))
window_stats$disparity_change_rate <- c(NA, diff(window_stats$morpho_disparity_mean) / diff(window_stats$avg_age_mean))

theme_evolution <- theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5), axis.title = element_text(size = 12),
        axis.text = element_text(size = 10), legend.title = element_text(size = 11),
        legend.text = element_text(size = 10), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.text = element_text(size = 11, face = "bold"))

p1 <- ggplot(morpho_data, aes(x = avg_age, y = morpho_disparity)) +
  geom_point(aes(color = distance_to_center), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, color = "red", linewidth = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  scale_color_viridis_c(name = "Distance to\nCenter") +
  labs(title = "Temporal Evolution of Morphological Disparity",
       subtitle = paste("GAM R² =", round(gam_summary$gam_r2[gam_summary$variable == "morpho_disparity"], 3),
                        ", Linear R² =", round(gam_summary$linear_r2[gam_summary$variable == "morpho_disparity"], 3)),
       x = "Time (Ma)", y = "Morphological Disparity") + scale_x_reverse() + theme_evolution

p2 <- ggplot(morpho_data, aes(x = avg_age, y = PCO1)) +
  geom_point(aes(color = specialization), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, color = "red", linewidth = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "dashed") +
  scale_color_viridis_c(name = "Specialization\nDegree") +
  labs(title = paste("PCO1 Temporal Evolution (", round(variance_explained[1], 1), "% variance)"),
       subtitle = paste("Kendall τ =", round(trend_results$kendall_tau[1], 3), ", p =", round(trend_results$kendall_p[1], 3)),
       x = "Time (Ma)", y = "PCO1") + scale_x_reverse() + theme_evolution

morpho_data$evolution_pattern <- "Stable"
for (i in 1:nrow(morpho_data)) {
  disparity <- morpho_data$morpho_disparity[i]
  rate <- ifelse(is.na(morpho_data$overall_rate[i]), 0, morpho_data$overall_rate[i])
  specialization <- morpho_data$specialization[i]
  high_disparity <- disparity > quantile(morpho_data$morpho_disparity, 0.75)
  high_rate <- rate > quantile(morpho_data$overall_rate, 0.75, na.rm = TRUE)
  high_special <- specialization > quantile(morpho_data$specialization, 0.75)
  if (high_disparity && high_rate) {
    morpho_data$evolution_pattern[i] <- "Rapid Radiation"
  } else if (high_special && !high_rate) {
    morpho_data$evolution_pattern[i] <- "Specialized Evolution"
  } else if (high_rate && !high_disparity) {
    morpho_data$evolution_pattern[i] <- "Directional Evolution"
  }
}

write.csv(morpho_data, "morphospace_evolution_analysis.csv", row.names = FALSE)
write.csv(trend_results, "evolution_trends_statistics.csv", row.names = FALSE)
