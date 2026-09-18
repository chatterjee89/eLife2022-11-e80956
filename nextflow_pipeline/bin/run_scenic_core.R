#!/usr/bin/env Rscript
# Core R-SCENIC run (gene filtering -> GENIE3 co-expression -> regulons via
# RcisTarget -> AUCell cell scoring), generalized from
# BPTF_mammary_GRN_analyses/02_SCENIC/01_run_SCENIC.R. Organism, cisTarget
# database paths, and filtering thresholds are parameterized instead of
# hardcoded to mm10/mgi.
# AI usage note: this pipeline step was written by Claude (Anthropic) and
# reviewed/validated by the author, generalizing a script the author had
# written for one fixed organism/dataset.

suppressPackageStartupMessages({
  library(optparse)
  library(loomR)
  library(SCENIC)
})

opt_list <- list(
  make_option("--loom", type = "character", help = "Expression matrix as a .loom file (genes x cells)"),
  make_option("--org", type = "character", default = "mgi", help = "SCENIC org code, e.g. mgi (mouse), hgnc (human), dmel (fly) [default %default]"),
  make_option("--db_dir", type = "character", help = "Directory containing the cisTarget .feather databases"),
  make_option("--db_10kb", type = "character", help = "cisTarget database filename: 10kb up/down of TSS"),
  make_option("--db_500bp", type = "character", help = "cisTarget database filename: 500bp up / 100bp down of TSS"),
  make_option("--dataset_title", type = "character", default = "SCENIC run"),
  make_option("--min_counts_per_cell", type = "double", default = 3, help = "Per-gene min counts threshold multiplier (see geneFiltering) [default %default]"),
  make_option("--min_pct_cells", type = "double", default = 0.01, help = "Minimum fraction of cells a gene must be detected in [default %default]"),
  make_option("--n_cores", type = "integer", default = 4),
  make_option("--seed", type = "integer", default = 123),
  make_option("--outdir", type = "character", default = "scenic")
)
opt <- parse_args(OptionParser(option_list = opt_list))

int_dir <- file.path(opt$outdir, "int")
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)
old_wd <- setwd(opt$outdir)
on.exit(setwd(old_wd), add = TRUE)

# ---------------------------------------------------------------------------
# 1. LOAD EXPRESSION MATRIX FROM LOOM
# ---------------------------------------------------------------------------
loom <- open_loom(file.path("..", opt$loom))
exprMat <- get_dgem(loom)
close_loom(loom)
saveRDS(exprMat, file.path("int", "exprMat.Rds"))

# ---------------------------------------------------------------------------
# 2. INITIALISE SCENIC OPTIONS
# ---------------------------------------------------------------------------
defaultDbNames[[opt$org]][1] <- opt$db_10kb
defaultDbNames[[opt$org]][2] <- opt$db_500bp

scenicOptions <- initializeScenic(
  org          = opt$org,
  dbDir        = opt$db_dir,
  dbs          = defaultDbNames[[opt$org]],
  datasetTitle = opt$dataset_title
)
scenicOptions@settings$nCores <- opt$n_cores
saveRDS(scenicOptions, file = file.path("int", "scenicOptions.Rds"))

# ---------------------------------------------------------------------------
# 3. GENE FILTERING -- retain genes expressed in >= min_pct_cells of cells
# with >= min_counts_per_cell counts in those cells.
# ---------------------------------------------------------------------------
genesKept <- geneFiltering(
  exprMat,
  scenicOptions    = scenicOptions,
  minCountsPerGene = opt$min_counts_per_cell * opt$min_pct_cells * ncol(exprMat),
  minSamples       = ncol(exprMat) * opt$min_pct_cells
)
exprMat_filtered <- exprMat[genesKept, ]
saveRDS(exprMat_filtered, file.path("int", "exprMat_filtered.Rds"))

# ---------------------------------------------------------------------------
# 4. GENIE3 CO-EXPRESSION NETWORK -- TF -> target-gene importance matrix by
# random forest regression. Log2-transforms expression first (standard
# SCENIC practice).
# ---------------------------------------------------------------------------
exprMat_filtered_log <- log2(exprMat_filtered + 1)
runCorrelation(exprMat_filtered_log, scenicOptions)
runGenie3(exprMat_filtered_log, scenicOptions, resumePreviousRun = FALSE)

# ---------------------------------------------------------------------------
# 5. STEP 1 -- co-expression modules from GENIE3 importance scores
# ---------------------------------------------------------------------------
scenicOptions <- readRDS(file.path("int", "scenicOptions.Rds"))
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- opt$seed
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

# ---------------------------------------------------------------------------
# 6. STEP 2 -- regulons via RcisTarget motif enrichment (prunes modules to
# genes with a TF-binding motif in their promoter/enhancer region)
# ---------------------------------------------------------------------------
scenicOptions <- readRDS(file.path("int", "scenicOptions.Rds"))
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- opt$seed
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

# ---------------------------------------------------------------------------
# 7. STEP 3 -- score cells (AUCell) for each regulon
# ---------------------------------------------------------------------------
scenicOptions <- readRDS(file.path("int", "scenicOptions.Rds"))
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- opt$seed

exprMat_log <- log2(exprMat + 1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file = file.path("int", "scenicOptions.Rds"))

n_regulons <- tryCatch(length(loadInt(scenicOptions, "regulons")), error = function(e) NA)
writeLines(c(
  "# plot_type: 'html'",
  "# section_name: 'R-SCENIC core run'",
  sprintf("<p>Organism: <b>%s</b>, cells: <b>%d</b>, genes kept after filtering: <b>%d</b> / %d</p>",
          opt$org, ncol(exprMat), nrow(exprMat_filtered), nrow(exprMat)),
  sprintf("<p>Regulons identified: <b>%s</b></p>", ifelse(is.na(n_regulons), "n/a", n_regulons))
), file.path("scenic_core_summary_mqc.txt"))

cat("SCENIC core run done.\n")
