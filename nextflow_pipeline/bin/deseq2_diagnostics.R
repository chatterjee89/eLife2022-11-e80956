#!/usr/bin/env Rscript
# Pre-processing / data-distribution diagnostics via DESeq2, generalized from
# the eLife2022-11-e80956 bulkSeq_code.R DESeq2 section. This step is
# diagnostics-only (normalization QC plots, PCA, a significance-based
# clustering heatmap) -- it is NOT the differential expression test the
# pipeline reports; that's edger_de.R.
# AI usage note: ported/generalized with Claude (Anthropic); reviewed by the
# author.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(gplots)
  library(RColorBrewer)
  library(knitr)
})

opt_list <- list(
  make_option("--counts", type = "character"),
  make_option("--sample_info", type = "character"),
  make_option("--sample_id_col", type = "character", default = "sample"),
  make_option("--group_columns", type = "character", default = "group",
              help = "Comma-separated sample_info columns combined (pasted with '.') into the DESeq2 design group [default %default]"),
  make_option("--pca_intgroup", type = "character", default = NULL,
              help = "Column used to color the PCA plot; defaults to the first group_columns entry"),
  make_option("--alpha", type = "double", default = 0.0001,
              help = "FDR cutoff for genes kept in the diagnostic clustering heatmap [default %default]"),
  make_option("--epsilon", type = "double", default = 1, help = "Pseudocount added before log2 in QC plots [default %default]"),
  make_option("--outdir", type = "character", default = "deseq2_diagnostics")
)
opt <- parse_args(OptionParser(option_list = opt_list))
epsilon <- opt$epsilon

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "svg"), showWarnings = FALSE)

# Load featureCounts output + sample metadata and build a DESeq2 dataset from
# the requested group columns. This pass exists to drive the diagnostic
# plots below -- it is not the DE test used downstream.
group_cols <- strsplit(opt$group_columns, ",")[[1]]
pca_intgroup <- if (is.null(opt$pca_intgroup)) group_cols[1] else opt$pca_intgroup

targets <- read.csv(opt$sample_info)
group <- factor(apply(targets[group_cols], 1, paste, collapse = "."))
t <- cbind(targets, group = group)

x <- read.csv(opt$counts, row.names = "Geneid")
x <- x[, as.character(t[[opt$sample_id_col]])]
count.table <- x

dds0 <- DESeqDataSetFromMatrix(countData = x, colData = t, design = ~group)
dds.norm <- estimateSizeFactors(dds0)

# Fixed color per group, reused across every diagnostic plot below.
palette <- colorRampPalette(brewer.pal(8, "Set2"))(nlevels(group))
col.sample <- setNames(palette, levels(group))
t$color <- col.sample[as.vector(t$group)]

# Raw vs. size-factor-normalized counts, as boxplots and density curves --
# confirms normalization brought the sample distributions into line before
# trusting any downstream comparison.
svg(file.path(opt$outdir, "svg", "qc_normalization.svg"), width = 10, height = 8)
par(mfrow = c(2, 2), cex.lab = 0.7)
boxplot(log2(counts(dds.norm) + epsilon), col = t$color, cex.axis = 0.7,
        las = 1, xlab = "log2(counts)", horizontal = TRUE, main = "Raw counts")
boxplot(log2(counts(dds.norm, normalized = TRUE) + epsilon), col = t$color, cex.axis = 0.7,
        las = 1, xlab = "log2(normalized counts)", horizontal = TRUE, main = "Normalized counts")
plotDensity(log2(counts(dds.norm) + epsilon), col = t$color,
            xlab = "log2(counts)", main = "Density: raw counts", cex.lab = 0.7, panel.first = grid())
plotDensity(log2(counts(dds.norm, normalized = TRUE) + epsilon), col = t$color,
            xlab = "log2(normalized counts)", main = "Density: normalized counts", cex.lab = 0.7, panel.first = grid())
dev.off()

# Per-sample summary stats -- a numeric sanity check alongside the plots above.
norm.counts <- counts(dds.norm, normalized = TRUE)
norm.counts.stats <- data.frame(
  min = apply(norm.counts, 2, min), mean = apply(norm.counts, 2, mean),
  median = apply(norm.counts, 2, median), max = apply(norm.counts, 2, max),
  zeros = apply(norm.counts == 0, 2, sum),
  percent.zeros = 100 * apply(norm.counts == 0, 2, sum) / nrow(norm.counts)
)
write.csv(norm.counts.stats, file.path(opt$outdir, "normalized_count_stats.csv"))

# Mean-variance relationship across genes -- RNA-seq counts are overdispersed
# relative to Poisson (variance > mean), which is why DESeq2/edgeR both model
# counts with a negative binomial instead.
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)
svg(file.path(opt$outdir, "svg", "mean_variance.svg"), width = 6, height = 6)
mean.var.col <- densCols(x = log2(mean.counts), y = log2(variance.counts))
plot(x = log2(mean.counts), y = log2(variance.counts), pch = 16, cex = 0.5, col = mean.var.col,
     main = "Mean-variance relationship", xlab = "Mean log2(normalized counts)",
     ylab = "Variance of log2(normalized counts)", panel.first = grid())
abline(a = 0, b = 1, col = "brown")
dev.off()

# Dispersion + Wald test used only to pick genes for the clustering heatmap
# below; the manuscript-equivalent reported DE gene list comes from edgeR
# (edger_de.R), not from this result.
dds.disp <- estimateDispersions(dds.norm)
svg(file.path(opt$outdir, "svg", "dispersion_estimates.svg"), width = 6, height = 6)
plotDispEsts(dds.disp)
dev.off()

wald.test <- nbinomWaldTest(dds.disp)
res.DESeq2 <- results(wald.test, alpha = opt$alpha, pAdjustMethod = "BH")

gene.kept <- rownames(res.DESeq2)[res.DESeq2$padj <= opt$alpha & !is.na(res.DESeq2$padj)]
count.table.kept <- log2(count.table + epsilon)[gene.kept, , drop = FALSE]

if (nrow(count.table.kept) >= 2) {
  svg(file.path(opt$outdir, "svg", "diagnostic_clustering_heatmap.svg"), width = 8, height = 10)
  heatmap.2(as.matrix(count.table.kept), scale = "row",
            hclust = function(x) hclust(x, method = "average"),
            distfun = function(x) as.dist((1 - cor(t(x))) / 2),
            trace = "none", density = "none", labRow = "", cexCol = 0.7)
  dev.off()
}

# Re-run the full DESeq2 pipeline in one call (rather than the separate
# estimateSizeFactors/estimateDispersions/nbinomWaldTest steps above) purely
# to get a variance-stabilized matrix for the PCA plot.
dds <- DESeq(dds0, betaPrior = FALSE)
vst <- varianceStabilizingTransformation(dds)
svg(file.path(opt$outdir, "svg", "pca.svg"), width = 6, height = 5)
print(plotPCA(vst, intgroup = pca_intgroup))
dev.off()

summary_txt <- file.path(opt$outdir, "deseq2_diagnostics_summary_mqc.txt")
writeLines(c(
  "# plot_type: 'html'",
  "# section_name: 'DESeq2 diagnostics (QC only, not the reported DE test)'",
  sprintf("<p>Samples: <b>%d</b>, groups: <b>%s</b></p>", ncol(x), paste(levels(group), collapse = ", ")),
  sprintf("<p>Genes at FDR &le; %.2g used for the clustering heatmap: <b>%d</b></p>", opt$alpha, length(gene.kept))
), summary_txt)

cat("DESeq2 diagnostics done.\n")
