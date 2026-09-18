#!/usr/bin/env Rscript
# Export a Seurat object carrying velocyto assays (RNA/spliced/unspliced/
# ambiguous, plus an "integrated" assay from clustering) to loom and h5ad for
# downstream RNA velocity analysis in scVelo (see scvelo_analysis.py).
# Generalized from the eLife2022-11-e80956 seurat_code.R "SCVELO PREPARATION"
# section. In the original, this step reads a separately prepared,
# already-subset velocity object (not derived within the same script from
# the ctrl/treatment integration above) -- this port keeps that same
# decoupling: --velocity_rds is any Seurat object with those assays, however
# it was produced (e.g. a cluster-of-interest subset of seurat_process.R's
# integrated_final.rds).
# AI usage note: ported/generalized with Claude (Anthropic); reviewed by the
# author.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(SeuratDisk)
})

opt_list <- list(
  make_option("--velocity_rds", type = "character", help = "Seurat .rds with RNA/spliced/unspliced/ambiguous(/integrated) assays"),
  make_option("--loom_out", type = "character", default = "velocity.loom"),
  make_option("--h5ad_out", type = "character", default = "velocity.h5ad"),
  make_option("--assays", type = "character", default = "RNA,spliced,unspliced,ambiguous,integrated",
              help = "Comma-separated assays to scale in turn before export; assays absent from the object are skipped [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

velgl <- readRDS(opt$velocity_rds)

# Scale each requested assay in turn so all are populated regardless of
# which one the exported file is read back into downstream.
for (a in strsplit(opt$assays, ",")[[1]]) {
  if (a %in% Assays(velgl)) {
    DefaultAssay(velgl) <- a
    velgl <- FindVariableFeatures(velgl)
    velgl <- ScaleData(velgl)
  } else {
    message(sprintf("Assay '%s' not present, skipping.", a))
  }
}
DefaultAssay(velgl) <- "RNA"

# Export to loom (scVelo's native format) ...
velgl.loom <- as.loom(velgl, filename = opt$loom_out, verbose = TRUE, overwrite = TRUE)
velgl.loom$close_all()

# ... and separately to h5ad via an intermediate h5Seurat file -- this is the
# file scvelo_analysis.py actually reads.
h5seurat_path <- sub("\\.h5ad$", ".h5Seurat", opt$h5ad_out)
SaveH5Seurat(velgl, filename = h5seurat_path, overwrite = TRUE)
Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
if (file.exists(sub("\\.h5Seurat$", ".h5ad", h5seurat_path)) &&
    sub("\\.h5Seurat$", ".h5ad", h5seurat_path) != opt$h5ad_out) {
  file.rename(sub("\\.h5Seurat$", ".h5ad", h5seurat_path), opt$h5ad_out)
}

cat("Velocity export done:", opt$loom_out, "and", opt$h5ad_out, "\n")
