# Single-cell transcriptomics identifies Keap1-Nrf2 regulated collective invasion in a *Drosophila* tumor model

**eLife 2022;11:e80956 · DOI: [10.7554/eLife.80956](https://doi.org/10.7554/eLife.80956)**

---

## Overview

This repository contains all analysis code for the study, which uses *Drosophila melanogaster* fat body as a tumor model. Tumors are induced by RNAi-mediated knockdown of *lgl* (lethal giant larvae) using the temperature-sensitive driver *tjTS*, and compared to GFP-expressing controls at two timepoints (24 h and 96 h post-induction). The study combines:

- **Whole-tissue bulk RNA-seq** to profile transcriptional changes at the tissue level
- **Single-cell RNA-seq (scRNA-seq)** to resolve cell-type-specific responses within the tumor
- **RNA velocity** to infer directional transcriptional dynamics in the invasive cell population

---

## Experimental design

| Sample | Genotype | Timepoint | Abbreviation |
|--------|----------|-----------|--------------|
| Control replicate 1 & 2 | *tjTS>GFP* (w1118) | 24 h | WT_24 |
| Control replicate 1 & 2 | *tjTS>GFP* (w1118) | 96 h | WT_96 |
| Tumor replicate 1 & 2 | *tjTS>lglRNAi* | 24 h | LglIR_24 |
| Tumor replicate 1 & 2 | *tjTS>lglRNAi* | 96 h | LglIR_96 |

scRNA-seq data are available on SRA: control (SRX7814226 / SRP321151), lglRNAi (LglKD_FC).

---

## Repository structure

```
├── preprocessing_code.sh   # Shell: STAR alignment, featureCounts, Cell Ranger, Velocyto
├── bulkSeq_code.R          # R: bulk RNA-seq QC, DESeq2 diagnostics, edgeR DGE
├── seurat_code.R           # R: scRNA-seq clustering and CCA integration (Seurat v4)
└── scvelo_code.py          # Python: RNA velocity analysis (scVelo)
```

---

## Analysis workflow

### Step 1 — Raw data preprocessing (`preprocessing_code.sh`)

**Bulk RNA-seq**
1. Build STAR genome index from *D. melanogaster* BDGP6.28 genome + custom GAL4/EGFP sequences
2. Align each of the 8 samples (40 threads, unique mappers only: `--outFilterMultimapNmax 1`)
3. Count reads per gene with `featureCounts` → `featurecounts_bulkSeq.csv`

**scRNA-seq**
1. Quantify spliced counts with **Cell Ranger 3.0.1** (`SC3Pv2` chemistry) for both control and lglRNAi samples
2. Generate spliced/unspliced count matrices with **Velocyto** for downstream RNA velocity analysis

---

### Step 2 — Bulk RNA-seq analysis (`bulkSeq_code.R`)

**QC and diagnostics (DESeq2)**
- Size-factor normalisation and count distribution visualisation (box plots, density plots)
- Mean–variance relationship and negative binomial distribution verification
- Dispersion estimation and `plotDispEsts`
- Hierarchical clustering heatmap of DEGs (padj ≤ 0.0001, Pearson-correlation distance, average linkage) → **Fig. 1e**
- PCA on variance-stabilised counts → **Fig. 1d**

**Differential gene expression (edgeR)**
- TMM normalisation → GLM quasi-likelihood F-test
- Key contrast: **96 h *lglRNAi* vs all other conditions** (the dominant axis of PC1)
- Effect-size threshold: `glmTreat` with |log2FC| ≥ log2(1.2)
- Output: `DGEgenes_LglRNAi96HvsOthers(PC1).csv` → **Supplementary Data 1**

---

### Step 3 — scRNA-seq clustering and integration (`seurat_code.R`)

**Per-sample processing (Seurat v4)**

| QC parameter | Control (w1118) | lglRNAi |
|-------------|-----------------|---------|
| Max UMI / spliced counts | 12,000 | 15,000 |
| Feature range | 650–1800 | 600–1800 |
| Max mitochondrial % | 10% | 4.5% |

Both samples undergo:
- Cell-cycle scoring (G2M and S scores) → regression of `CC.Difference`
- Variable feature selection (VST, 2,000 features)
- Scaling (regressing out `CC.Difference`, `percent.mt`, `nCount_RNA`)
- PCA (100 PCs) → clustering and UMAP with 60 PCs

**Data integration (CCA)**
- Canonical correlation analysis via `FindIntegrationAnchors` / `IntegrateData` (dims 1–50)
- Two-round integration: initial merge → re-split by sample → re-integrate on RNA assay
- Final clustering: resolution 1, UMAP (`n.neighbors = 100`, `min.dist = 0.50`)
- Removal of low-quality / doublet clusters after visual inspection

**scVelo preparation**
- Export the lglRNAi cluster-7 subset (invasive population) to `.loom` and `.h5ad` formats using `SeuratDisk`

---

### Step 4 — RNA velocity (`scvelo_code.py`)

Applied to **cluster 7 cells** (invasive population) exported from Seurat:

1. **Steady-state model** — initial velocity estimate and gene filtering (`min_r2 = 1`)
2. **Dynamical model** — kinetic parameter recovery via `recover_dynamics` for higher accuracy
3. Latent time estimation → pseudotime-ordered heatmap of top velocity genes
4. Velocity confidence scoring and gene ranking per cluster

---

## Required input files

| File | Used by | Description |
|------|---------|-------------|
| `~/cellranger_output/` | `preprocessing_code.sh` | Cell Ranger FASTQ input |
| `~/STAR_index/` | `preprocessing_code.sh` | Pre-built STAR genome index |
| `~/Drosophila_melanogaster.BDGP6.28.dna.toplevel_GAL4_EGFP.fa` | `preprocessing_code.sh` | Custom genome FASTA |
| `~/Drosophila_melanogaster.BDGP6.28.100.chr_filtered.GAL4.EGFP.gtf` | `preprocessing_code.sh` | Custom annotation GTF |
| `~/dm6_rmsk.gtf` | `preprocessing_code.sh` | RepeatMasker GTF for Velocyto |
| `featurecounts_bulkSeq.csv` | `bulkSeq_code.R` | featureCounts output |
| `~/sample_info.csv` | `bulkSeq_code.R` | Sample metadata with `Genotype`, `Time`, `PC1` columns |
| `~/cell_cycle_genes.txt` | `bulkSeq_code.R`, `seurat_code.R` | Cell cycle gene list (lines 1–68: G2/M; 69–124: S phase) |
| `~/Ctrl_FC/` | `seurat_code.R` | Cell Ranger output for control |
| `~/LglIR.loom` | `seurat_code.R` | Velocyto loom for lglRNAi sample |
| `lgl.h5ad` | `scvelo_code.py` | h5ad exported from Seurat (cluster 7 subset) |

---

## Dependencies

**R packages**
```r
# Bulk RNA-seq
library(DESeq2); library(edgeR); library(limma)
library(pheatmap); library(ggplot2); library(gplots); library(RColorBrewer)

# scRNA-seq
library(Seurat)      # v4
library(SeuratDisk)
```

**Python packages**
```python
import scvelo as scv   # RNA velocity
```

**Command-line tools**
- STAR 2.7
- Cell Ranger 3.0.1
- Velocyto
- Subread 2.0.0 (featureCounts)

---

## Citation

> Chatterjee D, Bhavna R, Bhatt A, Bhatt P, Bhatt S, et al.  
> *Single-cell transcriptomics identifies Keap1-Nrf2 regulated collective invasion in a Drosophila tumor model.*  
> **eLife** 2022;11:e80956. https://doi.org/10.7554/eLife.80956
