# AI usage note: this analysis code was written by the author. Claude Code
# (Anthropic) was used afterward, retrospectively, to organize the
# repository (file layout, README) — not to write or modify the analysis
# logic below. See AI_DISCLOSURE.md for the full disclosure.
#CONTROL: w1118
ctrl.data <- Read10X(data.dir = "~/Ctrl_FC/new_reference/")
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "Ctrl_FC", min.cells = 3, min.features = 200)
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^mt:")
ctrl[["percent.rp"]] <- PercentageFeatureSet(ctrl, pattern = "^Rp")
ctrl <- subset(ctrl, subset = nCount_RNA < 12000 & nFeature_RNA < 1800 & nFeature_RNA > 650 & percent.mt < 10)
ctrl <- NormalizeData(ctrl)
cc.genes <- readLines(con = "~/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
ctrl <- CellCycleScoring(ctrl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
ctrl$CC.Difference <- ctrl$S.Score - ctrl$G2M.Score
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
ctrl <- ScaleData(ctrl, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
ctrl <- RunPCA(ctrl, features = VariableFeatures(ctrl), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
ctrl <- FindNeighbors(ctrl, dims = 1:60)
ctrl <- FindClusters(ctrl, resolution = 1)
ctrl <- RunUMAP(ctrl, dims = 1:60, n.neighbors = 40, min.dist = 0.25)
DimPlot(ctrl, reduction = "umap", pt.size = 1.2, label = TRUE)
ctrl <- subset(ctrl, idents = c(19,17,8,13,15,9), invert = TRUE)

#SAMPLE: lglIR
ldat <- ReadVelocity(file = "~/LglIR.loom")
lgl <- as.Seurat(x = ldat)
lgl[["percent.mt"]] <- PercentageFeatureSet(lgl, pattern = "^mt:")
lgl[["percent.rp"]] <- PercentageFeatureSet(lgl, pattern = "^Rp")
lgl <- subset(lgl, subset = nCount_spliced < 15000 & nFeature_spliced < 1800 & nFeature_spliced > 600 & percent.mt < 4.5)
lgl <- NormalizeData(lgl)
cc.genes <- readLines(con = "~/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
lgl <- CellCycleScoring(lgl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
lgl$CC.Difference <- lgl$S.Score - lgl$G2M.Score
lgl <- FindVariableFeatures(lgl, selection.method = "vst", nfeatures = 2000)
lgl <- ScaleData(lgl, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_spliced"), verbose = TRUE)
lgl <- RunPCA(lgl, features = VariableFeatures(pbmc), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
lgl <- FindNeighbors(lgl, dims = 1:60)
lgl <- FindClusters(lgl, resolution = 1)
table(lgl@active.ident,lgl@meta.data$orig.ident)
lgl <- RunUMAP(lgl, dims = 1:60, n.neighbors = 40, min.dist = 0.20)
DimPlot(lgl, reduction = "umap", pt.size = 1.2, label = TRUE)

lgl <- subset(lgl, idents = c(11,8,16,12), invert = T)
lgl@assays$spliced -> lgl@assays$RNA
lgl@meta.data$nFeature_spliced -> lgl@meta.data$nFeature_RNA
lgl@meta.data$nCount_spliced -> lgl@meta.data$nCount_RNA
lgl@meta.data$orig.ident <- "LglIR_FC"

#DATA INTEGRATION
all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, lgl), anchor.features = 2000, dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
all.combined -> pbmc

all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
all.combined <- FindNeighbors(all.combined, dims = 1:60)
all.combined <- FindClusters(all.combined, resolution = 1.5)
table(all.combined@active.ident,all.combined@meta.data$orig.ident)
all.combined <- RunUMAP(all.combined, dims = 1:60, n.neighbors = 40, min.dist = 0.20)
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = F, group.by = "orig.ident")
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = T)

all.combined <- subset(all.combined, idents = c(24), invert = T)

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

all.combined <- RunPCA(all.combined, features = VariableFeatures(all.combined), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
all.combined <- FindNeighbors(all.combined, dims = 1:60)
all.combined <- FindClusters(all.combined, resolution = 1)
table(all.combined@active.ident,all.combined@meta.data$orig.ident)
all.combined <- RunUMAP(all.combined, dims = 1:60, n.neighbors = 100, min.dist = 0.50)
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = TRUE)
DimPlot(all.combined, reduction = "umap", pt.size = 1.2, label = F, group.by = "orig.ident")

#SCVELO PREPARATION
library(SeuratDisk)
velgl <- readRDS("~/lglIR_velocity.rds")

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
velgl.loom <- as.loom(velgl, filename = "~/cluster7ONLY_vel.loom", verbose = TRUE, overwrite = TRUE)
velgl.loom

SaveH5Seurat(velgl, filename = "lgl.h5Seurat", overwrite = TRUE)
Convert("lgl.h5Seurat", dest = "h5ad", overwrite = TRUE)

velgl.loom$close_all()
