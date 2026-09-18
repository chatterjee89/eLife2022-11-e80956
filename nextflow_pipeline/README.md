# nf-multiomics-de-grn

> AI tools were used in developing this repository — see [../AI_DISCLOSURE.md](../AI_DISCLOSURE.md) for details.

This is the Nextflow version of the analysis in the parent repo. If you're
here to redo **this repo's** analysis (the *Drosophila* lgl-RNAi tumor
model), you want `--mode bulk_de` and `--mode scrna_velocity` — see the
"Reproducing this repo's analysis via Nextflow" section in
[../README.md](../README.md) for exact commands with this repo's own file
names filled in. `--mode scenic`/`full` port a separate, unrelated pipeline
(mouse BPTF GRN inference) and don't apply to this repo's data; they're
included here because this Nextflow pipeline is a single combined tool
rather than being duplicated across two repos.

A generalized Nextflow (DSL2) pipeline combining three analyses that started
as separate, dataset-specific scripts:

| Mode              | Ported from                                                                 | What it does |
|-------------------|------------------------------------------------------------------------------|--------------|
| `bulk_de`         | `eLife2022-11-e80956/bulkSeq_code.R`                                        | DESeq2 QC/diagnostic plots (optional) + the edgeR differential expression test |
| `scrna_velocity`  | `eLife2022-11-e80956/seurat_code.R` + `scvelo_code.py`                      | Seurat QC/clustering/integration for a control-vs-treatment scRNA-seq design, plus an optional RNA velocity leg via scVelo |
| `scenic`          | `BPTF_mammary_GRN_analyses/02_SCENIC/{01_run_SCENIC.R,02_SCENIC_downstream_analysis.R}` | R-SCENIC gene regulatory network inference (gene filtering → GENIE3 → regulons → AUCell), with optional RSS/CSI downstream analysis |
| `full`            | all of the above                                                             | Runs all three subworkflows independently in the same execution |

These three started from unrelated studies (a *Drosophila* tumor model and a
mouse BPTF knockout), so `full` mode does not assume any biological
relationship between the inputs — it just runs all three subworkflows
side by side against whatever inputs you give each one.

## Important: pre-processed inputs only

This pipeline does **not** run alignment (STAR) or velocyto — those are
assumed to have already been run externally, matching what the original
scripts assumed. It starts from:

- `bulk_de`: a featureCounts-style count matrix + a sample metadata CSV
- `scrna_velocity`: a 10x Genomics matrix directory (control) + a velocyto
  `.loom` file (treatment); the RNA-velocity leg additionally needs a
  pre-assembled Seurat object carrying `RNA`/`spliced`/`unspliced`/
  `ambiguous` assays (this mirrors the original `seurat_code.R`, where the
  velocity-prep step reads its own separately prepared object rather than
  one built earlier in the same script)
- `scenic`: an expression matrix as a `.loom` file, plus cisTarget motif
  databases for your organism

## Usage

```bash
# Validate wiring with tiny synthetic inputs, no real compute:
nextflow run main.nf -profile test -stub-run --mode full

# Bulk RNA-seq DE only:
nextflow run main.nf -profile conda --mode bulk_de \
  --bulk_counts featurecounts.csv --bulk_sample_info sample_info.csv \
  --bulk_contrast_column PC1 --bulk_contrast_levels yes,no

# scRNA-seq + integration (velocity leg skipped unless --scrna_velocity_rds is given):
nextflow run main.nf -profile conda --mode scrna_velocity \
  --scrna_ctrl_10x_dir ctrl_10x/ --scrna_treatment_loom treatment.loom \
  --scrna_cell_cycle_genes cell_cycle_genes.txt

# SCENIC GRN inference:
nextflow run main.nf -profile conda --mode scenic \
  --scenic_loom expr.loom --scenic_org mgi \
  --scenic_cistarget_dir /path/to/cisTarget_databases \
  --scenic_db_10kb mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
  --scenic_db_500bp mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather

# All three, on an SGE cluster (mirrors the original 01_run_SCENIC.sh qsub setup):
nextflow run main.nf -profile conda,sge --mode full ...
```

Every threshold that was a hardcoded literal in the original scripts (QC
cutoffs, clustering resolutions, UMAP parameters, FDR/fold-change cutoffs,
cisTarget database paths, etc.) is now a `--param`; see `nextflow.config`
for the full list and its defaults, which match the original scripts'
values.

## Structure

```
main.nf                    # mode dispatch + subworkflow wiring
subworkflows/               # one workflow per mode (bulk_de, scrna_velocity, scenic)
modules/                    # one process per pipeline step, each a thin wrapper around bin/
bin/                        # the actual R/Python analysis logic (parameterized, callable standalone)
envs/                       # one conda environment per process group
conf/test.config            # tiny synthetic inputs + fast settings for -stub-run
assets/test_data/           # the synthetic inputs referenced by conf/test.config
```

## Environments

Each process declares its own conda environment (`envs/*.yml`), resolved
per-process rather than one shared environment — DESeq2/edgeR, Seurat,
scVelo, and R-SCENIC have overlapping-but-conflicting dependency trees, and
this avoids that conflict entirely. Run with `-profile conda` (or
`-profile conda,sge` on a cluster) to enable it.

## Reporting

Every step writes a small `*_mqc.txt`/plot pair in
[MultiQC custom-content](https://multiqc.info/docs/custom_content/) format
alongside its normal outputs; a final `MULTIQC` step (on by default, disable
with `--run_multiqc false`) aggregates whichever steps ran into one HTML
report under `<outdir>/report/`.

## Known deviation from the original scripts

The original `seurat_code.R` line building the treatment sample's PCA used
`VariableFeatures(pbmc)`, where `pbmc` was a variable not yet assigned at
that point in the script (fixed in the eLife repo's own commit history).
This port uses the corrected, self-consistent form throughout.
