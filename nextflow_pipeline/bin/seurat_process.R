#!/usr/bin/env Rscript
# scRNA-seq QC, clustering, and two-round integration for a 2-condition
# (control vs. treatment) design, generalized from the eLife2022-11-e80956
# seurat_code.R control/treatment/integration sections. Deliberately kept
# 2-condition (not N-sample) to match the structure of the original data;
# paths, thresholds, and clustering parameters are configurable.
#
# The original script's RunPCA(lgl, features = VariableFeatures(pbmc), ...)
# referenced a variable not yet assigned at that point in the script (fixed
# upstream in the eLife repo to VariableFeatures(lgl) -- see its commit
# history); this port uses the corrected, self-consistent form throughout.
#
# AI usage note: ported/generalized with Claude (Anthropic) from code the
# author originally wrote for one specific dataset; reviewed by the author.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(velocyto.R)
})

opt_list <- list(
  make_option("--ctrl_10x_dir", type = "character", help = "10x Genomics matrix directory for the control sample"),
  make_option("--ctrl_sample_name", type = "character", default = "Ctrl"),
  make_option("--treatment_loom", type = "character", help = "velocyto .loom file for the treatment sample"),
  make_option("--treatment_sample_name", type = "character", default = "Treatment"),

  make_option("--cell_cycle_genes", type = "character", help = "Newline-delimited gene list: first n_g2m_genes are G2M, next n_s_genes are S"),
  make_option("--n_g2m_genes", type = "integer", default = 68),
  make_option("--n_s_genes", type = "integer", default = 56),

  make_option("--mito_pattern", type = "character", default = "^mt:"),
  make_option("--ribo_pattern", type = "character", default = "^Rp"),

  make_option("--ctrl_max_count", type = "double", default = 12000),
  make_option("--ctrl_max_feature", type = "double", default = 1800),
  make_option("--ctrl_min_feature", type = "double", default = 650),
  make_option("--ctrl_max_mito", type = "double", default = 10),

  make_option("--treat_max_count", type = "double", default = 15000),
  make_option("--treat_max_feature", type = "double", default = 1800),
  make_option("--treat_min_feature", type = "double", default = 600),
  make_option("--treat_max_mito", type = "double", default = 4.5),

  make_option("--npcs", type = "integer", default = 100),
  make_option("--dims", type = "integer", default = 60, help = "Number of PCs used for neighbors/UMAP throughout"),
  make_option("--nfeatures", type = "integer", default = 2000, help = "Number of HVGs per FindVariableFeatures call"),

  make_option("--ctrl_resolution", type = "double", default = 1),
  make_option("--ctrl_umap_neighbors", type = "integer", default = 40),
  make_option("--ctrl_umap_min_dist", type = "double", default = 0.25),
  make_option("--ctrl_drop_clusters", type = "character", default = "", help = "Comma-separated cluster ids to drop after inspecting the ctrl UMAP"),

  make_option("--treat_resolution", type = "double", default = 1),
  make_option("--treat_umap_neighbors", type = "integer", default = 40),
  make_option("--treat_umap_min_dist", type = "double", default = 0.20),
  make_option("--treat_drop_clusters", type = "character", default = "", help = "Comma-separated cluster ids to drop after inspecting the treatment UMAP"),

  make_option("--anchor_features", type = "integer", default = 2000),
  make_option("--anchor_dims", type = "integer", default = 50),

  make_option("--integration1_resolution", type = "double", default = 1.5),
  make_option("--integration1_umap_neighbors", type = "integer", default = 40),
  make_option("--integration1_umap_min_dist", type = "double", default = 0.20),
  make_option("--integration1_drop_clusters", type = "character", default = "", help = "Comma-separated cluster ids to drop after the first integration pass"),

  make_option("--integration2_resolution", type = "double", default = 1),
  make_option("--integration2_umap_neighbors", type = "integer", default = 100),
  make_option("--integration2_umap_min_dist", type = "double", default = 0.50),

  make_option("--outdir", type = "character", default = "seurat_process")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "svg"), showWarnings = FALSE)
svg_plot <- function(name, expr, width = 7, height = 6) {
  svg(file.path(opt$outdir, "svg", paste0(name, ".svg")), width = width, height = height)
  print(expr)
  dev.off()
}
drop_ids <- function(csv) {
  csv <- trimws(csv)
  if (csv == "") return(integer(0))
  as.integer(strsplit(csv, ",")[[1]])
}

cc.genes <- readLines(opt$cell_cycle_genes)
g2m.genes <- cc.genes[seq_len(opt$n_g2m_genes)]
s.genes <- cc.genes[(opt$n_g2m_genes + 1):(opt$n_g2m_genes + opt$n_s_genes)]

cluster_and_embed <- function(obj, resolution, umap_neighbors, umap_min_dist, count_col) {
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = opt$nfeatures)
  obj <- ScaleData(obj, vars.to.regress = c("CC.Difference", "percent.mt", count_col), verbose = TRUE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = opt$npcs, nfeature.print = 10, ndims.print = 1:5, verbose = TRUE)
  obj <- FindNeighbors(obj, dims = 1:opt$dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = 1:opt$dims, n.neighbors = umap_neighbors, min.dist = umap_min_dist)
  obj
}

# ---------------------------------------------------------------------------
# CONTROL SAMPLE
# ---------------------------------------------------------------------------
ctrl.data <- Read10X(data.dir = opt$ctrl_10x_dir)
ctrl <- CreateSeuratObject(counts = ctrl.data, project = opt$ctrl_sample_name, min.cells = 3, min.features = 200)
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = opt$mito_pattern)
ctrl[["percent.rp"]] <- PercentageFeatureSet(ctrl, pattern = opt$ribo_pattern)
ctrl <- subset(ctrl, subset = nCount_RNA < opt$ctrl_max_count & nFeature_RNA < opt$ctrl_max_feature &
                 nFeature_RNA > opt$ctrl_min_feature & percent.mt < opt$ctrl_max_mito)
ctrl <- NormalizeData(ctrl)
ctrl <- CellCycleScoring(ctrl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ctrl$CC.Difference <- ctrl$S.Score - ctrl$G2M.Score
ctrl <- cluster_and_embed(ctrl, opt$ctrl_resolution, opt$ctrl_umap_neighbors, opt$ctrl_umap_min_dist, "nCount_RNA")
svg_plot("ctrl_umap", DimPlot(ctrl, reduction = "umap", pt.size = 1.2, label = TRUE))
ctrl <- subset(ctrl, idents = drop_ids(opt$ctrl_drop_clusters), invert = TRUE)

# ---------------------------------------------------------------------------
# TREATMENT SAMPLE (velocyto loom -> spliced counts repointed to RNA assay)
# ---------------------------------------------------------------------------
ldat <- ReadVelocity(file = opt$treatment_loom)
treat <- as.Seurat(x = ldat)
treat[["percent.mt"]] <- PercentageFeatureSet(treat, pattern = opt$mito_pattern)
treat[["percent.rp"]] <- PercentageFeatureSet(treat, pattern = opt$ribo_pattern)
treat <- subset(treat, subset = nCount_spliced < opt$treat_max_count & nFeature_spliced < opt$treat_max_feature &
                  nFeature_spliced > opt$treat_min_feature & percent.mt < opt$treat_max_mito)
treat <- NormalizeData(treat)
treat <- CellCycleScoring(treat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
treat$CC.Difference <- treat$S.Score - treat$G2M.Score
treat <- cluster_and_embed(treat, opt$treat_resolution, opt$treat_umap_neighbors, opt$treat_umap_min_dist, "nCount_spliced")

write.csv(as.data.frame.matrix(table(treat@active.ident, treat@meta.data$orig.ident)),
          file.path(opt$outdir, "treatment_cluster_by_sample.csv"))
svg_plot("treatment_umap", DimPlot(treat, reduction = "umap", pt.size = 1.2, label = TRUE))
treat <- subset(treat, idents = drop_ids(opt$treat_drop_clusters), invert = TRUE)

# Repoint the default RNA slots at the spliced counts so `treat`'s assay
# structure matches `ctrl`'s ahead of integration (IntegrateData expects a
# standard RNA assay, not the velocyto spliced/unspliced/ambiguous split).
treat@assays$spliced -> treat@assays$RNA
treat@meta.data$nFeature_spliced -> treat@meta.data$nFeature_RNA
treat@meta.data$nCount_spliced -> treat@meta.data$nCount_RNA
treat@meta.data$orig.ident <- opt$treatment_sample_name

# ---------------------------------------------------------------------------
# INTEGRATION -- pass 1 (catch any remaining low-quality cluster)
# ---------------------------------------------------------------------------
all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, treat), anchor.features = opt$anchor_features, dims = 1:opt$anchor_dims)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:opt$anchor_dims)
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)

all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = opt$npcs, nfeature.print = 10, ndims.print = 1:5, verbose = TRUE)
all.combined <- FindNeighbors(all.combined, dims = 1:opt$dims)
all.combined <- FindClusters(all.combined, resolution = opt$integration1_resolution)
write.csv(as.data.frame.matrix(table(all.combined@active.ident, all.combined@meta.data$orig.ident)),
          file.path(opt$outdir, "integration1_cluster_by_sample.csv"))
all.combined <- RunUMAP(all.combined, dims = 1:opt$dims, n.neighbors = opt$integration1_umap_neighbors, min.dist = opt$integration1_umap_min_dist)
svg_plot("integration1_umap_by_sample", DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = FALSE, group.by = "orig.ident"))
svg_plot("integration1_umap_clusters", DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = TRUE))

all.combined <- subset(all.combined, idents = drop_ids(opt$integration1_drop_clusters), invert = TRUE)

# ---------------------------------------------------------------------------
# INTEGRATION -- pass 2 (re-split, re-normalize fresh, re-integrate on
# cleaned data; this is the final clustering used downstream)
# ---------------------------------------------------------------------------
split.combined <- SplitObject(all.combined, split.by = "orig.ident")
ctrl <- split.combined[[opt$ctrl_sample_name]]
treat <- split.combined[[opt$treatment_sample_name]]
DefaultAssay(ctrl) <- "RNA"
DefaultAssay(treat) <- "RNA"
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = opt$nfeatures)
treat <- FindVariableFeatures(treat, selection.method = "vst", nfeatures = opt$nfeatures)
ctrl <- NormalizeData(ctrl)
treat <- NormalizeData(treat)

all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, treat), anchor.features = opt$anchor_features, dims = 1:opt$anchor_dims)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:opt$anchor_dims)
DefaultAssay(all.combined) <- "integrated"
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)

all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = opt$npcs, nfeature.print = 10, ndims.print = 1:5, verbose = TRUE)
all.combined <- FindNeighbors(all.combined, dims = 1:opt$dims)
all.combined <- FindClusters(all.combined, resolution = opt$integration2_resolution)
write.csv(as.data.frame.matrix(table(all.combined@active.ident, all.combined@meta.data$orig.ident)),
          file.path(opt$outdir, "integration2_cluster_by_sample.csv"))
all.combined <- RunUMAP(all.combined, dims = 1:opt$dims, n.neighbors = opt$integration2_umap_neighbors, min.dist = opt$integration2_umap_min_dist)
svg_plot("final_umap_clusters", DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = TRUE))
svg_plot("final_umap_by_sample", DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = FALSE, group.by = "orig.ident"))

saveRDS(all.combined, file.path(opt$outdir, "integrated_final.rds"))

summary_txt <- file.path(opt$outdir, "seurat_process_summary_mqc.txt")
writeLines(c(
  "# plot_type: 'html'",
  "# section_name: 'Seurat QC / clustering / integration'",
  sprintf("<p>Final integrated object: <b>%d</b> cells, <b>%d</b> clusters (resolution %.2f)</p>",
          ncol(all.combined), nlevels(Idents(all.combined)), opt$integration2_resolution)
), summary_txt)

cat("Seurat processing done:", ncol(all.combined), "cells in final integrated object.\n")
