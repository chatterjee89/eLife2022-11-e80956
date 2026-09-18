# AI usage note: this analysis code was written by the author. Claude Code
# (Anthropic) was used afterward, retrospectively, to organize the
# repository (file layout, README) — not to write or modify the analysis
# logic below. See AI_DISCLOSURE.md for the full disclosure.
# A generalized, parameterized Nextflow version of this script (`--mode
# scrna_velocity`) is available in nextflow_pipeline/ -- see its README, or
# the "Reproducing this analysis via Nextflow" section of the repo-root
# README.
#CONTROL: w1118
# Load w1118 control 10x Genomics counts and build the base Seurat object.
ctrl.data <- Read10X(data.dir = "~/Ctrl_FC/new_reference/")
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "Ctrl_FC", min.cells = 3, min.features = 200)

# QC metrics: mitochondrial and ribosomal read fraction (Drosophila gene
# prefixes "mt:" / "Rp"), used to filter dying/low-quality cells below.
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^mt:")
ctrl[["percent.rp"]] <- PercentageFeatureSet(ctrl, pattern = "^Rp")

# Drop low-count, high-mito, and likely-doublet (high-feature) cells.
ctrl <- subset(ctrl, subset = nCount_RNA < 12000 & nFeature_RNA < 1800 & nFeature_RNA > 650 & percent.mt < 10)
ctrl <- NormalizeData(ctrl)

# Cell-cycle scoring so its effect can be regressed out during scaling --
# without this, clusters can form around cycling state rather than cell type.
cc.genes <- readLines(con = "~/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
ctrl <- CellCycleScoring(ctrl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ctrl$CC.Difference <- ctrl$S.Score - ctrl$G2M.Score

# Standard Seurat clustering: HVGs -> scale (regressing out cell cycle,
# mito%, and depth) -> PCA -> graph-based clustering -> UMAP.
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
ctrl <- ScaleData(ctrl, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
ctrl <- RunPCA(ctrl, features = VariableFeatures(ctrl), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
ctrl <- FindNeighbors(ctrl, dims = 1:60)
ctrl <- FindClusters(ctrl, resolution = 1)
ctrl <- RunUMAP(ctrl, dims = 1:60, n.neighbors = 40, min.dist = 0.25)
DimPlot(ctrl, reduction = "umap", pt.size = 1.2, label = TRUE)

# Remove clusters identified as low-quality/non-epithelial/doublet after
# visual inspection of the UMAP above.
ctrl <- subset(ctrl, idents = c(19,17,8,13,15,9), invert = TRUE)

#SAMPLE: lglIR
# Load the lgl-RNAi tumor sample from its velocyto loom file (carries
# spliced/unspliced/ambiguous counts needed later for RNA velocity) and wrap
# it as a Seurat object.
ldat <- ReadVelocity(file = "~/LglIR.loom")
lgl <- as.Seurat(x = ldat)
lgl[["percent.mt"]] <- PercentageFeatureSet(lgl, pattern = "^mt:")
lgl[["percent.rp"]] <- PercentageFeatureSet(lgl, pattern = "^Rp")

# Same QC-filter / normalize / cell-cycle-score pipeline as the control
# sample, but against the spliced-count slots (this object came from
# velocyto) and with thresholds tuned separately for this sample.
lgl <- subset(lgl, subset = nCount_spliced < 15000 & nFeature_spliced < 1800 & nFeature_spliced > 600 & percent.mt < 4.5)
lgl <- NormalizeData(lgl)
cc.genes <- readLines(con = "~/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
lgl <- CellCycleScoring(lgl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
lgl$CC.Difference <- lgl$S.Score - lgl$G2M.Score
lgl <- FindVariableFeatures(lgl, selection.method = "vst", nfeatures = 2000)
lgl <- ScaleData(lgl, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_spliced"), verbose = TRUE)
lgl <- RunPCA(lgl, features = VariableFeatures(lgl), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
lgl <- FindNeighbors(lgl, dims = 1:60)
lgl <- FindClusters(lgl, resolution = 1)

# Cluster-by-sample cross-tab -- quick check of cluster composition before
# deciding what to drop.
table(lgl@active.ident,lgl@meta.data$orig.ident)
lgl <- RunUMAP(lgl, dims = 1:60, n.neighbors = 40, min.dist = 0.20)
DimPlot(lgl, reduction = "umap", pt.size = 1.2, label = TRUE)

# Same low-quality/doublet cluster removal as for ctrl, based on this UMAP.
lgl <- subset(lgl, idents = c(11,8,16,12), invert = T)

# Repoint the default RNA slots at the spliced counts so `lgl`'s assay
# structure matches `ctrl`'s ahead of integration (IntegrateData expects a
# standard RNA assay, not the velocyto spliced/unspliced/ambiguous split).
lgl@assays$spliced -> lgl@assays$RNA
lgl@meta.data$nFeature_spliced -> lgl@meta.data$nFeature_RNA
lgl@meta.data$nCount_spliced -> lgl@meta.data$nCount_RNA
lgl@meta.data$orig.ident <- "LglIR_FC"

#DATA INTEGRATION
# Find shared anchor cells between ctrl and lgl and use them to integrate the
# two samples into a common space (removes sample-of-origin/batch effects
# ahead of joint clustering).
all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, lgl), anchor.features = 2000, dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
all.combined -> pbmc

# First-pass joint clustering, run mainly to spot any remaining
# low-quality/doublet cluster that survived per-sample QC above.
all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
all.combined <- FindNeighbors(all.combined, dims = 1:60)
all.combined <- FindClusters(all.combined, resolution = 1.5)
table(all.combined@active.ident,all.combined@meta.data$orig.ident)
all.combined <- RunUMAP(all.combined, dims = 1:60, n.neighbors = 40, min.dist = 0.20)
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = F, group.by = "orig.ident")
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = T)

# Drop the one remaining low-quality cluster found in this first-pass
# integration.
all.combined <- subset(all.combined, idents = c(24), invert = T)

# Split the once-integrated object back into per-sample pieces, then redo
# normalization/HVG selection on each (now QC'd) sample fresh, and
# re-integrate -- a cleaner second integration pass now that only
# good-quality cells remain on both sides.
split.combined <- SplitObject(all.combined, split.by = "orig.ident")
split.combined
ctrl <- split.combined$Ctrl_FC
lgl <- split.combined$LglIR_FC
DefaultAssay(ctrl) <- "RNA"
DefaultAssay(lgl) <- "RNA"
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
lgl <- FindVariableFeatures(lgl, selection.method = "vst", nfeatures = 2000)
ctrl <- NormalizeData(ctrl)
lgl <- NormalizeData(lgl)
all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, lgl), anchor.features = 2000, dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)
DefaultAssay(all.combined) <- "integrated"
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)

# Final joint clustering used for the paper's UMAP/cluster calls: same
# PCA -> neighbors -> clusters -> UMAP steps as before, now on the cleaned,
# re-integrated object (slightly different resolution/UMAP neighbor
# parameters than the first pass).
all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
all.combined <- FindNeighbors(all.combined, dims = 1:60)
all.combined <- FindClusters(all.combined, resolution = 1)
table(all.combined@active.ident,all.combined@meta.data$orig.ident)
all.combined <- RunUMAP(all.combined, dims = 1:60, n.neighbors = 100, min.dist = 0.50)
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = TRUE)
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = F, group.by = "orig.ident")

#SCVELO PREPARATION
# Export the finalized cluster assignments/embedding, together with the
# spliced/unspliced/ambiguous counts from a separately saved velocity
# object, to loom/h5ad format for downstream RNA velocity analysis in
# Python (see scvelo_code.py).
library(SeuratDisk)
velgl <- readRDS("~/lglIR_velocity.rds")

# Scale each assay in turn (RNA, then the three velocyto count layers, then
# the integrated assay) so all are populated regardless of which one the
# exported file is read back into downstream.
DefaultAssay(velgl) <- "RNA"
velgl <- FindVariableFeatures(velgl)
velgl <- ScaleData(velgl)
DefaultAssay(velgl) <- "spliced"
velgl <- FindVariableFeatures(velgl)
velgl <- ScaleData(velgl)
DefaultAssay(velgl) <- "unspliced"
velgl <- FindVariableFeatures(velgl)
velgl <- ScaleData(velgl)
DefaultAssay(velgl) <- "ambiguous"
velgl <- FindVariableFeatures(velgl)
velgl <- ScaleData(velgl)
DefaultAssay(velgl) <- "integrated"
velgl <- FindVariableFeatures(velgl)
velgl <- ScaleData(velgl)
DefaultAssay(velgl) <- "RNA"

# Export to loom (scVelo's native format) ...
velgl.loom <- as.loom(velgl, filename = "~/cluster7ONLY_vel.loom", verbose = TRUE, overwrite = TRUE)
velgl.loom

# ... and separately to h5ad via an intermediate h5Seurat file -- this is
# the file scvelo_code.py actually reads in as "lgl.h5ad".
SaveH5Seurat(velgl, filename = "lgl.h5Seurat", overwrite = TRUE)
Convert("lgl.h5Seurat", dest = "h5ad", overwrite = TRUE)

velgl.loom$close_all()
