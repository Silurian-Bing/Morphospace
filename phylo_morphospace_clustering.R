req_files <- c("T29h tree.nex", "ragesave.R", "taxon_age_data.csv", "morphospace.xlsx")
if (!all(file.exists(req_files))) stop("Required data files missing")

library(ape);library(readxl);library(dplyr);library(RColorBrewer);library(ggplot2);library(gridExtra)
library(cluster);library(factoextra);library(ggtree);library(phytools);library(viridis)
setwd("F:/HB/DataBase/Morphospace/")

tr <- read.nexus("T29h tree.nex"); if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
tree_tx <- tr$tip.label; source("ragesave.R"); age_dt <- read.csv("taxon_age_data.csv")
xl_file <- "morphospace.xlsx"; sht_names <- excel_sheets(xl_file)

calc_phylo_metrics <- function(tree, data, method = "K") {
  if (require(phytools, quietly = TRUE)) { return(phylosig(tree, data, method = method)) }
  return(NA)
}

get_optimal_k <- function(data, max_k, min_size) {
  sil_scores <- numeric(max_k-1)
  for (k in 2:max_k) {
    km_res <- kmeans(data, centers = k, nstart = 50)
    grp_sizes <- table(km_res$cluster)
    if (min(grp_sizes) >= min_size) {
      dt_dist <- dist(data); sil <- silhouette(km_res$cluster, dt_dist)
      sil_scores[k-1] <- mean(sil[, 3])
    } else { sil_scores[k-1] <- -1 }
  }
  valid_k <- which(sil_scores > 0) + 1
  return(if(length(valid_k) == 0) 3 else valid_k[which.max(sil_scores[valid_k-1])])
}

all_morph <- list(); all_tx <- c()
for (sht in sht_names) {
  curr_dt <- read_excel(xl_file, sheet = sht)
  if (ncol(curr_dt) > 1) {
    tx_col <- curr_dt[[1]]; val_rows <- !is.na(tx_col) & tx_col != ""
    if (sum(val_rows) > 0) {
      sht_dt <- curr_dt[val_rows, ]; names(sht_dt)[1] <- "taxon"
      sht_dt$geo_stage <- sht; all_morph[[sht]] <- sht_dt; all_tx <- c(all_tx, sht_dt$taxon)
    }
  }
}

common_tx <- intersect(intersect(tree_tx, age_dt$taxon), unique(all_tx))
if (length(common_tx) < 15) stop("Insufficient taxa for analysis")

if (length(all_morph) > 0) {
  all_colnms <- lapply(all_morph, names); comm_cols <- Reduce(intersect, all_colnms)
  final_morph <- data.frame()
  for (sht_dt in all_morph) { sht_sub <- sht_dt[, comm_cols]; final_morph <- rbind(final_morph, sht_sub) }
  final_morph <- final_morph[final_morph$taxon %in% common_tx, ]
}

morph_feat <- final_morph[, !names(final_morph) %in% c("taxon", "geo_stage")]
morph_mtx <- matrix(NA, nrow = nrow(morph_feat), ncol = ncol(morph_feat))
colnames(morph_mtx) <- names(morph_feat); rownames(morph_mtx) <- final_morph$taxon

val_cols <- logical(ncol(morph_feat))
for (i in 1:ncol(morph_feat)) {
  col_dt <- suppressWarnings(as.numeric(morph_feat[[i]]))
  if (sum(!is.na(col_dt)) >= 5 && var(col_dt, na.rm = TRUE) > 0) {
    morph_mtx[, i] <- col_dt; val_cols[i] <- TRUE
  }
}

morph_mtx <- morph_mtx[, val_cols, drop = FALSE]
val_rows <- rowSums(!is.na(morph_mtx)) >= 3; morph_mtx <- morph_mtx[val_rows, , drop = FALSE]
avail_tx <- rownames(morph_mtx); final_tx <- intersect(intersect(avail_tx, age_dt$taxon), tree_tx)
final_morph_cln <- morph_mtx[final_tx, , drop = FALSE]

for (i in 1:ncol(final_morph_cln)) {
  na_idx <- is.na(final_morph_cln[, i])
  if (any(na_idx)) { col_mn <- mean(final_morph_cln[, i], na.rm = TRUE); final_morph_cln[na_idx, i] <- col_mn }
}

final_age <- age_dt[age_dt$taxon %in% final_tx, ]
final_age <- final_age[match(final_tx, final_age$taxon), ]; final_age$avg_age <- (final_age$FAD + final_age$LAD) / 2
final_tr <- keep.tip(tr, final_tx)

pca_res <- prcomp(final_morph_cln, scale. = TRUE, center = TRUE)
var_exp <- summary(pca_res)$importance[2, ] * 100; cumvar <- cumsum(var_exp)
n_pcs_80 <- which(cumvar >= 80)[1]; n_pcs <- min(6, ncol(pca_res$x)); pc_scr <- pca_res$x[, 1:n_pcs]

pc_dt <- data.frame(taxon = rownames(pc_scr), PC1 = pc_scr[, 1], PC2 = if(n_pcs > 1) pc_scr[, 2] else 0,
                    PC3 = if(n_pcs > 2) pc_scr[, 3] else 0, avg_age = final_age$avg_age,
                    FAD = final_age$FAD, LAD = final_age$LAD, stringsAsFactors = FALSE)

n_pcs_clust <- min(n_pcs_80, n_pcs); pca_cols_clust <- paste0("PC", 1:n_pcs_clust)
pca_for_clust <- pc_dt[, pca_cols_clust, drop = FALSE]; pca_scaled <- scale(pca_for_clust)

set.seed(123); min_grp_size <- max(5, length(final_tx) %/% 8); max_k <- min(6, length(final_tx) %/% min_grp_size)
optimal_k <- get_optimal_k(pca_scaled, max_k, min_grp_size)

clust_res <- kmeans(pca_scaled, centers = optimal_k, nstart = 100)
pc_dt$morpho_grp <- paste0("Group_", clust_res$cluster)

pca_dist <- dist(pca_scaled); final_sil <- silhouette(clust_res$cluster, pca_dist)
avg_sil <- mean(final_sil[, 3])

grp_stats <- pc_dt %>% group_by(morpho_grp) %>%
  summarise(n_taxa = n(), pc1_mean = mean(PC1), pc1_sd = sd(PC1), pc2_mean = mean(PC2), pc2_sd = sd(PC2),
            age_mean = mean(avg_age), age_range = max(avg_age) - min(avg_age), .groups = 'drop')

pc1_time_cor <- cor.test(pc_dt$avg_age, pc_dt$PC1, method = "spearman")
pc2_time_cor <- cor.test(pc_dt$avg_age, pc_dt$PC2, method = "spearman")

time_bins <- cut(pc_dt$avg_age, breaks = 4, labels = c("Youngest", "Young", "Old", "Oldest"))
pc_dt$time_bin <- time_bins

disp_by_time <- pc_dt %>% group_by(time_bin) %>%
  summarise(n_taxa = n(), pc1_var = var(PC1), pc2_var = var(PC2), 
            total_var = var(PC1) + var(PC2), .groups = 'drop')

tree_order <- match(final_tr$tip.label, pc_dt$taxon)
phylo_sig_pc1 <- calc_phylo_metrics(final_tr, pc_dt$PC1[tree_order])
phylo_sig_pc2 <- calc_phylo_metrics(final_tr, pc_dt$PC2[tree_order])

grp_colors <- brewer.pal(min(max(3, optimal_k), 11), "Set2")
if (optimal_k > 11) grp_colors <- rainbow(optimal_k)
names(grp_colors) <- sort(unique(pc_dt$morpho_grp))

theme_evo <- theme_minimal() + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                                     axis.title = element_text(size = 12), panel.border = element_rect(color = "black", fill = NA))

p1 <- ggplot(pc_dt, aes(x = PC1, y = PC2)) + geom_point(aes(color = morpho_grp, size = avg_age), alpha = 0.8) +
  stat_ellipse(aes(color = morpho_grp), level = 0.68, linetype = "dashed", alpha = 0.7) +
  scale_color_manual(values = grp_colors, name = "Morpho Group") +
  scale_size_continuous(name = "Age (Ma)", range = c(2, 5)) +
  labs(title = "Morphospace Clustering Patterns", x = paste("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste("PC2 (", round(var_exp[2], 1), "%)")) + theme_evo

p2 <- ggplot(pc_dt, aes(x = PC1, y = PC2)) + geom_point(aes(color = avg_age), size = 3, alpha = 0.8) +
  scale_color_viridis(name = "Age\n(Ma)", option = "plasma") +
  labs(title = "Temporal Morphospace Trajectory", x = paste("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste("PC2 (", round(var_exp[2], 1), "%)")) + theme_evo

p3 <- ggplot(pc_dt, aes(x = avg_age, y = PC1)) + geom_point(aes(color = morpho_grp), size = 3, alpha = 0.8) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed", span = 0.8, alpha = 0.3) +
  scale_color_manual(values = grp_colors, name = "Morpho Group") +
  labs(title = "Temporal Evolution of PC1", x = "Average Age (Ma)", y = paste("PC1 (", round(var_exp[1], 1), "%)")) + theme_evo

disp_clean <- disp_by_time[!is.na(disp_by_time$total_var), ]
p4 <- ggplot(disp_clean, aes(x = time_bin, y = total_var)) + geom_col(fill = "steelblue", alpha = 0.7, width = 0.6) +
  geom_text(aes(label = round(total_var, 2)), vjust = -0.5, size = 4) +
  labs(title = "Morphological Disparity Through Time", x = "Time Period", y = "Disparity (PC1+PC2 variance)") +
  theme_evo + theme(axis.text.x = element_text(angle = 45, hjust = 1))

tree_tips <- final_tr$tip.label; pc_dt_tree <- pc_dt[match(tree_tips, pc_dt$taxon), ]
matched_tips <- !is.na(pc_dt_tree$taxon)

tree_dt <- data.frame(tip_label = tree_tips, morpho_grp = ifelse(matched_tips, pc_dt_tree$morpho_grp, "Unclassified"),
                      stringsAsFactors = FALSE)
extended_colors <- grp_colors; if (any(!matched_tips)) extended_colors <- c(grp_colors, "Unclassified" = "gray70")

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

if (require(ggtree, quietly = TRUE)) {
  tree_plot <- ggtree(final_tr) + geom_tiplab(size = 2, hjust = -0.1) +
    xlim(0, max(node.depth.edgelength(final_tr)) * 1.3) +
    geom_tippoint(color = extended_colors[tree_dt$morpho_grp], size = 3) +
    ggtitle("Phylogenetic Relationships and Morphological Groups") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  print(tree_plot)
}

write.csv(pc_dt, "phylo_morphospace_analysis.csv", row.names = FALSE)
write.csv(grp_stats, "morphological_groups_stats.csv", row.names = FALSE)