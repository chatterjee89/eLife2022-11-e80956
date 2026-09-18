#!/usr/bin/env Rscript
# SCENIC downstream analysis: AUC binarization, Regulon Specificity Score
# (RSS, Jensen-Shannon divergence), and Connection Specificity Index (CSI).
# Generalized from
# BPTF_mammary_GRN_analyses/02_SCENIC/02_SCENIC_downstream_analysis.R.
#
# Scope note: the original script also includes several illustrative,
# dataset-specific examples (a hardcoded 3-TF heatmap subset, a single-TF
# AUCell_plotTSNE call, a single hardcoded regulon's CSI ranking plot). Those
# are one-off exploratory calls rather than generalizable pipeline steps, so
# this port keeps the reusable computational core (binarization, RSS, CSI,
# CSI module activity) and their aggregate outputs/plots, and leaves the
# one-off examples out.
#
# AI usage note: this pipeline step was written by Claude (Anthropic) and
# reviewed/validated by the author, generalizing a script the author had
# written for one fixed dataset.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(philentropy)
  library(svMisc)
  library(SCENIC)
  library(Seurat)
  library(pheatmap)
  library(viridis)
  library(ggplot2)
  library(ggrepel)
  library(ggridges)
  library(ComplexHeatmap)
})

opt_list <- list(
  make_option("--scenic_options_rds", type = "character", help = "int/scenicOptions.Rds from run_scenic_core.R"),
  make_option("--seurat_rds", type = "character", help = "Seurat object whose metadata supplies the cell-type/cluster grouping"),
  make_option("--cell_type_column", type = "character", default = "seurat_clusters", help = "Seurat metadata column used as the RSS/CSI grouping variable [default %default]"),
  make_option("--rss_threshold", type = "double", default = 0.4, help = "RSS cutoff for the cell-type-specific regulon heatmap [default %default]"),
  make_option("--csi_nclust", type = "integer", default = 10, help = "Number of CSI modules to cut the clustering into [default %default]"),
  make_option("--outdir", type = "character", default = "scenic_downstream")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "svg"), showWarnings = FALSE)
svg_plot <- function(name, expr, width = 8, height = 8) {
  svg(file.path(opt$outdir, "svg", paste0(name, ".svg")), width = width, height = height)
  print(expr)
  dev.off()
}

# =============================================================================
# CUSTOM FUNCTIONS (ported as-is from the original, computation unchanged)
# =============================================================================

# Per-regulon AUC binarization threshold via k-means (k=2, faster and less
# conservative than AUCell_exploreThresholds). Threshold = max AUC of the
# low-activity cluster; cells with AUC=0 are excluded before clustering.
auc_thresh_kmeans <- function(regulonAUC) {
  regulons <- rownames(regulonAUC@assays@data@listData$AUC)
  kmeans_thresholds <- list()
  for (regulon_no in seq_along(regulons)) {
    progress(regulon_no)
    regulon <- regulons[regulon_no]
    df <- data.frame(
      auc = regulonAUC@assays@data@listData$AUC[regulon, ],
      cells = names(regulonAUC@assays@data@listData$AUC[regulon, ]),
      regulon = regulon
    )
    df <- df %>% subset(auc > 0)
    if (nrow(df) < 2) { kmeans_thresholds[[regulon]] <- 0; next }
    km <- kmeans(df$auc, centers = 2)
    df$cluster <- as.factor(km$cluster)
    cluster1_max <- max(subset(df, cluster == 1)$auc)
    cluster2_max <- max(subset(df, cluster == 2)$auc)
    if (cluster1_max > cluster2_max) {
      df <- df %>%
        mutate(cluster = gsub(2, 3, cluster)) %>%
        mutate(cluster = gsub(1, 2, cluster)) %>%
        mutate(cluster = gsub(3, 1, cluster))
    }
    df_sub <- df %>% arrange(desc(auc)) %>% subset(cluster == 1)
    kmeans_thresholds[[regulon]] <- df_sub[1, ]$auc
  }
  kmeans_thresholds
}

# Converts continuous AUC scores to binary (0/1) using per-regulon thresholds.
binarize_regulons <- function(regulonAUC, thresholds) {
  binary_regulon_list <- list()
  for (regulon_no in seq_along(names(thresholds))) {
    progress(regulon_no)
    regulon <- names(thresholds)[regulon_no]
    auc_df <- data.frame(
      auc = regulonAUC@assays@data@listData$AUC[regulon, ],
      cells = names(regulonAUC@assays@data@listData$AUC[regulon, ])
    )
    auc_df <- auc_df %>%
      mutate(regulon = if_else(auc >= thresholds[[regulon]], 1, 0)) %>%
      select(-auc)
    colnames(auc_df) <- c("cells", regulon)
    binary_regulon_list[[regulon]] <- auc_df
  }
  binary_regulon_list
}

# Regulon Specificity Score: 1 - sqrt(Jensen-Shannon divergence) between a
# regulon's binary activity distribution and a cell type's membership vector.
# Higher RSS -> more cell-type-specific regulon.
calculate_rss <- function(metadata, binary_regulons, cell_type_column) {
  cell_types <- unique(metadata[, cell_type_column])
  regulons <- rownames(binary_regulons)
  jsd_matrix_ct <- data.frame(regulon = character(), cell_type = character(), RSS = numeric())
  for (ct in unique(cell_types)) {
    for (regulon in regulons) {
      regulon_vec <- binary_regulons[regulon, ]
      regulon_vec_sum <- sum(regulon_vec)
      if (regulon_vec_sum > 0) {
        regulon_norm <- regulon_vec / regulon_vec_sum
        genotype_vec <- metadata[colnames(binary_regulons), ] %>%
          mutate(cell_class = if_else(get(cell_type_column) == ct, 1, 0))
        genotype_norm <- genotype_vec$cell_class / sum(genotype_vec$cell_class)
        dist_df <- rbind(regulon_norm, genotype_norm)
        jsd_divergence <- suppressMessages(philentropy::JSD(dist_df))
        rss <- 1 - sqrt(jsd_divergence)
        jsd_matrix_ct <- rbind(jsd_matrix_ct, data.frame(regulon = regulon, cell_type = ct, RSS = rss[1]))
      }
    }
  }
  jsd_matrix_ct %>% arrange(desc(RSS))
}

# Connection Specificity Index: for each regulon pair, the fraction of all
# other regulons with lower Pearson correlation to BOTH members of the pair
# than they have to each other. CSI=1 -> maximally co-specific pair.
calculate_csi <- function(regulonAUC, calc_extended = FALSE, verbose = FALSE) {
  compare_pcc <- function(vector_of_pcc, pcc) {
    pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
    if (pcc_larger == length(vector_of_pcc)) return(0)
    length(vector_of_pcc)
  }
  calc_csi <- function(reg, reg2, pearson_cor) {
    test_cor <- pearson_cor[reg, reg2]
    total_n <- ncol(pearson_cor)
    pcc_sub <- subset(pearson_cor, rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)
    sums <- apply(pcc_sub, MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
    length(sums[sums == nrow(pcc_sub)]) / total_n
  }
  regulonAUC_sub <- regulonAUC@assays@data@listData$AUC
  regulonAUC_sub <- if (calc_extended) {
    subset(regulonAUC_sub, grepl("extended", rownames(regulonAUC_sub)))
  } else {
    subset(regulonAUC_sub, !grepl("extended", rownames(regulonAUC_sub)))
  }
  pearson_cor <- cor(t(regulonAUC_sub))
  regulon_names <- rownames(pearson_cor)
  n <- length(regulon_names)
  csi_regulons <- data.frame(regulon_1 = character(n * n), regulon_2 = character(n * n), CSI = numeric(n * n))
  f <- 0
  for (reg in regulon_names) {
    if (verbose) print(reg)
    for (reg2 in regulon_names) {
      f <- f + 1
      csi_regulons[f, ] <- list(reg, reg2, calc_csi(reg, reg2, pearson_cor))
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  csi_regulons
}

plot_csi_modules <- function(csi_df, nclust = 10, font_size_regulons = 6) {
  csi_mat <- csi_df %>%
    spread(regulon_2, CSI) %>%
    tibble::column_to_rownames("regulon_1") %>%
    as.matrix()
  pheatmap(csi_mat, show_colnames = FALSE, color = viridis(n = 10),
           cutree_cols = nclust, cutree_rows = nclust, fontsize_row = font_size_regulons,
           cluster_cols = TRUE, cluster_rows = TRUE, treeheight_row = 300, treeheight_col = 300,
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
}

plot_rss_ranking <- function(rss_df, ggrepel_force = 1, ggrepel_point_padding = 0.2, top_genes = 4) {
  rss_df <- rss_df %>% subset(!grepl("extended", regulon))
  rss_sub <- rss_df %>% group_by(cell_type) %>% mutate(rank = order(order(RSS, decreasing = TRUE)))
  ggplot(rss_sub, aes(rank, RSS, label = regulon)) +
    geom_point(color = "grey20", size = 2) +
    geom_point(data = subset(rss_sub, rank < top_genes), color = "red", size = 2) +
    geom_text_repel(data = subset(rss_sub, rank < top_genes), force = ggrepel_force, point.padding = ggrepel_point_padding) +
    labs(x = "Rank", y = "RSS") +
    facet_wrap(~cell_type)
}

# =============================================================================
# LOAD DATA
# =============================================================================
scenicOptions <- readRDS(opt$scenic_options_rds)
seurat <- readRDS(opt$seurat_rds)
cellInfo <- data.frame(row.names = colnames(seurat))
cellInfo[[opt$cell_type_column]] <- seurat@meta.data[[opt$cell_type_column]]

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]

# =============================================================================
# AVERAGE REGULON ACTIVITY PER CLUSTER
# =============================================================================
regulonActivity_byCluster <- sapply(
  split(rownames(cellInfo), cellInfo[[opt$cell_type_column]]),
  function(cells) rowMeans(getAUC(regulonAUC)[, cells, drop = FALSE])
)
regulonActivity_byCluster_Scaled <- t(scale(t(regulonActivity_byCluster), center = TRUE, scale = TRUE))
write.csv(regulonActivity_byCluster, file.path(opt$outdir, "regulon_activity_by_cluster.csv"))
svg_plot("regulon_activity_heatmap", ComplexHeatmap::Heatmap(regulonActivity_byCluster_Scaled, name = "Regulon activity"))

# =============================================================================
# BINARIZATION
# =============================================================================
kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)
binary_regulons <- binarize_regulons(regulonAUC, kmeans_thresholds)
binary_matrix <- binary_regulons %>%
  purrr::reduce(left_join, by = "cells") %>%
  tibble::column_to_rownames("cells")
binary_matrix_t <- as.matrix(t(binary_matrix))
write.csv(binary_matrix_t, file.path(opt$outdir, "binary_regulon_matrix.csv"))

# =============================================================================
# REGULON SPECIFICITY SCORE (RSS)
# =============================================================================
rss_df <- calculate_rss(metadata = cellInfo, binary_regulons = binary_matrix_t, cell_type_column = opt$cell_type_column)
write.csv(rss_df, file.path(opt$outdir, "regulon_specificity_score.csv"), row.names = FALSE)
svg_plot("rss_ranking_by_cell_type", plot_rss_ranking(rss_df), width = 10, height = 8)

rss_nona <- subset(rss_df, RSS > 0)
svg_plot("rss_distribution_ridge",
         ggplot(rss_nona, aes(RSS, cell_type, fill = cell_type)) +
           geom_density_ridges(scale = 5, alpha = 0.75) +
           geom_vline(xintercept = 0.1) +
           theme(legend.position = "none"))

rss_wide <- rss_df %>% spread(cell_type, RSS)
rownames(rss_wide) <- rss_wide$regulon
rss_wide <- rss_wide[, 2:ncol(rss_wide), drop = FALSE]
rss_specific <- rss_wide[apply(rss_wide, 1, function(x) any(x > opt$rss_threshold, na.rm = TRUE)), , drop = FALSE]
if (nrow(rss_specific) > 0) {
  write.csv(rss_specific, file.path(opt$outdir, sprintf("rss_specific_regulons_gt_%.2f.csv", opt$rss_threshold)))
}

# =============================================================================
# CONNECTION SPECIFICITY INDEX (CSI)
# =============================================================================
regulons_csi <- calculate_csi(regulonAUC, calc_extended = FALSE)
write.csv(regulons_csi, file.path(opt$outdir, "connection_specificity_index.csv"), row.names = FALSE)
svg(file.path(opt$outdir, "svg", "csi_modules_heatmap.svg"), width = 12, height = 12)
plot_csi_modules(regulons_csi, nclust = opt$csi_nclust, font_size_regulons = 8)
dev.off()

writeLines(c(
  "# plot_type: 'html'",
  "# section_name: 'SCENIC downstream (RSS / CSI)'",
  sprintf("<p>Regulons: <b>%d</b>, cell types (%s): <b>%d</b></p>",
          nrow(regulonActivity_byCluster), opt$cell_type_column, ncol(regulonActivity_byCluster)),
  sprintf("<p>Regulons with RSS &gt; %.2f in at least one cell type: <b>%d</b></p>", opt$rss_threshold, nrow(rss_specific))
), file.path(opt$outdir, "scenic_downstream_summary_mqc.txt"))

cat("SCENIC downstream analysis done.\n")
