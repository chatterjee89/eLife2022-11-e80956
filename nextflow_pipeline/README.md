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

## Step-by-step guide for first-time users

No command-line or Nextflow experience assumed. If you're redoing this
repo's own analysis, the two commands you actually need are in the
[repo-root README](../README.md#reproducing-this-analysis-via-nextflow),
pre-filled with this repo's file names — this section explains what those
commands mean and how to get to the point of running them.

### 1. What you're installing (one-time setup, a few minutes)

Two pieces of software, neither of which is R, Python, Seurat, or DESeq2 —
the pipeline installs all of *those* for you automatically, per step, the
first time it needs them.

- **Nextflow** — the "conductor" that reads this pipeline and runs each
  step in the right order, keeping a log of what succeeded and what
  didn't. Install: <https://www.nextflow.io/docs/latest/install.html>
  (on Mac/Linux, it's usually one `curl` command).
- **conda** or **mamba** — a package manager that builds a self-contained
  toolbox of software for each step (so the DESeq2 step and the Seurat step
  can each get exactly the R packages they need, without conflicting with
  each other). If you already have Anaconda/Miniconda/Miniforge installed
  for other bioinformatics work, you already have this.

Check both are installed by opening a terminal and running:
```bash
nextflow -version
conda -version
```
If either prints "command not found," revisit its install step above
before continuing.

### 2. Where to run commands

Everything below is run from a **terminal** (Terminal.app on Mac, or your
cluster's SSH terminal), not by double-clicking any file. Navigate into the
folder where you downloaded/cloned this repository — the one that contains
this `nextflow_pipeline/` folder alongside `bulkSeq_code.R`, `seurat_code.R`,
etc. — for example:
```bash
cd ~/Downloads/eLife2022-11-e80956
```
Every command in this guide is run from *that* folder (the repo root, one
level above `nextflow_pipeline/`), not from inside `nextflow_pipeline/`
itself — that's why the commands below start with `nextflow run
nextflow_pipeline/main.nf` rather than just `nextflow run main.nf`.

### 3. Try it safely first — no real data, no long wait

Before pointing this at your actual sequencing data, run this once to
confirm Nextflow and the pipeline itself are wired up correctly. It uses
tiny placeholder files and skips all the real computation, so it finishes
in seconds rather than hours:
```bash
nextflow run nextflow_pipeline/main.nf -profile test -stub-run --mode full
```
If this prints `Pipeline completed successfully` at the end (skim past the
colorful progress lines above it), your setup is good and you can move on
to real data. If it errors out, that's an installation problem to fix now,
before you spend hours waiting on a real run only to hit the same error.

### 4. Run it on your real data

Make sure the input files this step needs are sitting in your repo-root
folder (same ones the original R/Python scripts needed — see "Required
input files" in the [repo-root README](../README.md)), then run the
relevant command from [there](../README.md#reproducing-this-analysis-via-nextflow).
The first real run will be slower than usual because conda is building
each step's toolbox for the first time — subsequent runs reuse those and
start immediately.

While it runs, Nextflow prints one line per step, updating in place, e.g.:
```
executor >  local (2)
[a1/2b3c4d] process > EDGER_DE (bulk_de)       [100%] 1 of 1 ✔
[f5/6g7h8i] process > DESEQ2_DIAGNOSTICS (bulk_de) [100%] 1 of 1 ✔
```
A ✔ means that step finished; a ✘ means it failed (see troubleshooting
below). This is normal to watch run for anywhere from minutes (bulk DE) to
hours (Seurat integration, SCENIC's GENIE3 step) depending on your data size
and machine.

### 5. Where your results end up

Everything lands in a new `results/` folder, created automatically inside
whichever directory you ran the command from — you don't need to create it
yourself, and nothing overwrites your original scripts or data. Inside,
outputs are grouped by mode and step, e.g. `results/bulk_de/edger_de/` holds
the differential expression gene table, `results/scrna_velocity/seurat_process/`
holds the clustering plots and the integrated object, and so on. If
`--run_multiqc` wasn't disabled (it's on by default), there's also a single
`results/report/multiqc_report.html` you can open in a browser for a
one-page summary of every step that ran.

### 6. If something goes wrong

- **Re-running after a fix**: add `-resume` to your command
  (`nextflow run nextflow_pipeline/main.nf -resume ...`). Nextflow will
  skip every step that already finished successfully and only re-run the
  one that failed (and anything downstream of it) — you don't lose
  completed work by fixing a typo and trying again.
- **Finding out *why* a step failed**: the terminal output names the failed
  step and gives a folder path like `work/a1/2b3c4d...` — inside that
  folder, `.command.log` (or `.command.err`) has the actual error message
  from the underlying R/Python script, same as if you'd run it by hand.
- **First run seems stuck on "creating conda environment"**: this is
  normal and can take several minutes per step the very first time; it's
  downloading and installing software, not hung. Only worry if it's been
  stuck with no disk/network activity for a long time.
- **Still stuck**: every parameter this guide's commands use (and many more
  for fine-tuning) is documented in `nextflow.config` in this folder, with
  the value each one defaults to.

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
