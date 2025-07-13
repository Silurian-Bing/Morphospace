req_file <- "morphospace.xlsx"
if (!file.exists(req_file)) stop("Required morphospace Excel file missing")

library(readxl);library(dplyr);library(tidyr);library(vegan);library(ggplot2);library(ggforce)
library(RColorBrewer);library(ggrepel);library(viridis);library(ggrastr)
setwd("F:/HB/DataBase/Morphospace/")

process_morpho_pca <- function(morph_data, meta_data) {
  if (nrow(morph_data) < 2 || ncol(morph_data) < 2) stop("Insufficient data for PCA")
  pca_res <- rda(morph_data, scale = TRUE)
  pca_scores <- scores(pca_res, display = "sites", choices = 1:min(5, ncol(morph_data)))
  return(list(result = pca_res, scores = as.data.frame(pca_scores)))
}

create_morpho_plots <- function(plot_data, pc1_var, pc2_var, color_pal, sheet_names) {
  p1 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = TimePeriod)) +
    geom_point(alpha = 0.8, size = 2.5) + scale_color_manual(values = color_pal, name = "Time Period") +
    labs(title = "Morphospace Evolution (PCA)", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)"),
         caption = "Points represent genera, colored by time period") +
    theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5, face = "bold")) + coord_fixed()
  
  p2 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = TimePeriod, fill = TimePeriod)) +
    geom_point(alpha = 0.6, size = 2) + ggforce::geom_mark_ellipse(aes(group = TimePeriod), alpha = 0.15) +
    scale_color_manual(values = color_pal) + scale_fill_manual(values = color_pal) +
    labs(title = "Morphospace Evolution with Ellipses", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)")) +
    theme_bw(base_size = 14) + guides(fill = "none") + coord_fixed()
  
  centroids <- plot_data %>% group_by(TimePeriod) %>%
    summarise(mean_PC1 = mean(PC1, na.rm = TRUE), mean_PC2 = mean(PC2, na.rm = TRUE), .groups = 'drop') %>%
    mutate(TimePeriod = factor(TimePeriod, levels = sheet_names, ordered = TRUE)) %>% arrange(TimePeriod)
  
  p3 <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = TimePeriod), alpha = 0.5, size = 2) +
    geom_path(data = centroids, aes(x = mean_PC1, y = mean_PC2), color = "black", linewidth = 1,
              arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
    geom_point(data = centroids, aes(x = mean_PC1, y = mean_PC2, fill = TimePeriod), size = 4, shape = 21, color = "black") +
    scale_color_manual(values = color_pal) + scale_fill_manual(values = color_pal) +
    labs(title = "Morphospace Evolution with Centroid Trajectory", x = paste0("PC1 (", pc1_var, "%)"), y = paste0("PC2 (", pc2_var, "%)")) +
    theme_bw(base_size = 14) + coord_fixed()
  
  return(list(points = p1, ellipses = p2, trajectory = p3))
}

file_path <- "morphospace.xlsx"
sht_names <- excel_sheets(file_path); time_periods <- factor(sht_names, levels = sht_names, ordered = TRUE)

all_dt_list <- lapply(sht_names, function(sht) {
  df <- tryCatch({ read_excel(path = file_path, sheet = sht, range = cell_cols("A:AV")) },
                 error = function(e) { warning(paste("Error reading sheet:", sht)); return(NULL) })
  
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return(NULL)
  if (ncol(df) != 48) warning(paste("Sheet:", sht, "- Read", ncol(df), "columns, expected 48"))
  
  df$TimePeriod <- sht; colnames(df)[1] <- "Genus"; df$Genus <- as.character(df$Genus)
  trait_cols <- names(df)[2:48]
  
  if (length(trait_cols) > 0) {
    df <- df %>% mutate(across(all_of(trait_cols), .fns = ~ {
      num_col <- suppressWarnings(as.numeric(as.character(.)))
      tidyr::replace_na(num_col, 0)
    }))
  }
  return(df)
})

all_dt_list <- all_dt_list[!sapply(all_dt_list, is.null)]
if (length(all_dt_list) == 0) stop("No data could be read successfully")

comb_dt <- bind_rows(all_dt_list)
if(nrow(comb_dt) == 0) stop("Combined data has 0 rows")

expected_trait_cols <- names(all_dt_list[[1]])[2:48]
char_cols <- intersect(expected_trait_cols, names(comb_dt))

morph_dt <- comb_dt %>% select(all_of(char_cols))
meta_dt <- comb_dt %>% select(Genus, TimePeriod)

if (any(is.na(morph_dt))) morph_dt[is.na(morph_dt)] <- 0

zero_var_cols <- which(apply(morph_dt, 2, var) == 0)
if (length(zero_var_cols) > 0) {
  morph_dt_complete <- morph_dt[, -zero_var_cols]
} else {
  morph_dt_complete <- morph_dt
}

if (nrow(meta_dt) != nrow(morph_dt_complete)) stop("Row mismatch between morpho data and metadata")

pca_res_list <- process_morpho_pca(morph_dt_complete, meta_dt)
pca_result <- pca_res_list$result; pca_scores_df <- pca_res_list$scores

plot_dt <- bind_cols(meta_dt, pca_scores_df)
plot_dt$TimePeriod <- factor(plot_dt$TimePeriod, levels = sht_names, ordered = TRUE)

pca_summary <- summary(pca_result)
exp_var <- pca_summary$cont$importance["Proportion Explained", ]
pc1_var <- round(exp_var["PC1"] * 100, 1); pc2_var <- round(exp_var["PC2"] * 100, 1)

n_periods <- length(unique(plot_dt$TimePeriod))
color_pal <- viridis_pal(option = "D")(n_periods)

plots <- create_morpho_plots(plot_dt, pc1_var, pc2_var, color_pal, sht_names)

ggsave("morphospace_pca_points.png", plots$points, width = 8, height = 6, dpi = 300)
ggsave("morphospace_pca_ellipses.png", plots$ellipses, width = 8, height = 6, dpi = 300)
ggsave("morphospace_pca_trajectory.png", plots$trajectory, width = 9, height = 7, dpi = 300)

pca_loadings <- scores(pca_result, display = "species", choices = 1:2)
if (nrow(pca_loadings) > 0) {
  top_pc1 <- head(pca_loadings[order(-abs(pca_loadings[, "PC1"])), , drop = FALSE])
  top_pc2 <- head(pca_loadings[order(-abs(pca_loadings[, "PC2"])), , drop = FALSE])
}

write.csv(plot_dt, "morphospace_pca_data.csv", row.names = FALSE)