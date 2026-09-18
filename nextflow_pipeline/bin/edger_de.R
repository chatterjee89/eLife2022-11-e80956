#!/usr/bin/env Rscript
# Differential expression via edgeR, generalized from the eLife2022-11-e80956
# bulkSeq_code.R edgeR section (CPM filter -> TMM -> GLM dispersion -> glmTreat).
# AI usage note: ported/generalized with Claude (Anthropic) from code the
# author originally wrote for a single fixed contrast; reviewed by the author.

suppressPackageStartupMessages({
  library(optparse)
  library(edgeR)
})

opt_list <- list(
  make_option("--counts", type = "character", help = "featureCounts-style count matrix CSV (genes x samples, Geneid rownames)"),
  make_option("--sample_info", type = "character", help = "Sample metadata CSV; must contain a sample id column and the contrast column"),
  make_option("--sample_id_col", type = "character", default = "sample", help = "Column in sample_info matching count matrix column names [default %default]"),
  make_option("--contrast_col", type = "character", help = "Column in sample_info defining the two groups to contrast (e.g. PC1)"),
  make_option("--contrast_levels", type = "character", help = "numerator,denominator levels of contrast_col, e.g. yes,no"),
  make_option("--min_cpm", type = "double", default = 1.0, help = "Minimum CPM to keep a gene [default %default]"),
  make_option("--min_samples", type = "integer", default = 2, help = "Minimum number of samples passing min_cpm to keep a gene [default %default]"),
  make_option("--lfc_threshold", type = "double", default = 1.2, help = "Fold-change threshold tested by glmTreat (linear scale, e.g. 1.2 = 1.2-fold) [default %default]"),
  make_option("--fdr", type = "double", default = 0.05, help = "BH-adjusted p-value cutoff for the reported gene table [default %default]"),
  make_option("--outdir", type = "character", default = "edger_de")
)
opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

targets <- read.csv(opt$sample_info)
levels_split <- strsplit(opt$contrast_levels, ",")[[1]]
numerator <- trimws(levels_split[1])
denominator <- trimws(levels_split[2])

group <- factor(targets[[opt$contrast_col]], levels = c(denominator, numerator))
t <- cbind(targets, group = group)

x <- read.csv(opt$counts, row.names = "Geneid")
# Align count matrix columns to sample_info row order via the sample id column.
x <- x[, as.character(t[[opt$sample_id_col]])]

genetable <- data.frame(gene_id = rownames(x))
y <- DGEList(counts = x, group = group, genes = genetable)

# Filter to genes with CPM > min_cpm in at least min_samples samples before
# normalizing -- standard edgeR practice, removes genes too lowly expressed
# to estimate dispersion reliably.
cpm_mat <- cpm(y)
countcheck <- cpm_mat > opt$min_cpm
keep <- which(rowSums(countcheck) >= opt$min_samples)
y <- y[keep, ]

y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)

# GLM dispersion estimation (common -> trended -> tagwise) and a
# quasi-likelihood fit -- edgeR's standard route to a robust per-gene
# dispersion estimate with a small sample size.
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmQLFit(y, design)

contrast_formula <- setNames(list(paste0(numerator, "-", denominator)), "contrast")
my.contrasts <- makeContrasts(contrasts = contrast_formula$contrast, levels = design)

# glmTreat tests against a fold-change threshold rather than just
# "significantly different from zero".
qlf <- glmTreat(fit, contrast = my.contrasts[, 1], lfc = log2(opt$lfc_threshold))

de_calls <- decideTestsDGE(qlf, adjust.method = "BH", p.value = opt$fdr)
tab <- topTags(qlf, n = Inf, sort.by = "PValue", p.value = opt$fdr)

out_csv <- file.path(opt$outdir, sprintf("DGEgenes_%s_vs_%s.csv", numerator, denominator))
write.csv(tab, file = out_csv)

summary_txt <- file.path(opt$outdir, "edger_de_summary_mqc.txt")
writeLines(c(
  "# plot_type: 'html'",
  "# section_name: 'edgeR differential expression'",
  sprintf("<p>Contrast: <b>%s vs %s</b> (column: %s)</p>", numerator, denominator, opt$contrast_col),
  sprintf("<p>Genes tested after filtering: <b>%d</b> / %d</p>", nrow(y), nrow(x)),
  sprintf("<p>Significant genes (FDR &lt; %.3g, |FC| &ge; %.2fx): <b>%d</b></p>", opt$fdr, opt$lfc_threshold, nrow(tab$table))
), summary_txt)

cat(sprintf("edgeR DE done: %d significant genes written to %s\n", nrow(tab$table), out_csv))
