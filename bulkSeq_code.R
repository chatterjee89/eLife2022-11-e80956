# AI usage note: this analysis code was written by the author. Claude Code
# (Anthropic) was used afterward, retrospectively, to organize the
# repository (file layout, README) — not to write or modify the analysis
# logic below. See AI_DISCLOSURE.md for the full disclosure.
library(edgeR)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(reticulate)
library(ggplot2)
library(airway)
library(affy)
library(magrittr)
library(gplots)
library(DESeq2)
library(lazyeval)
library(knitr)
library(gplots)

#DESEQ2 FOR PRE-PROCESSING AND DATA DISTRIBUTION DIAGNOSIS
targets<-read.csv("~/sample_info.csv")
group<-factor(paste(targets$Genotype,targets$Time,sep="."))
t <- cbind(targets,group=group)
x <- read.csv("~/featurecounts_bulkSeq.csv", row.names="Geneid")
x -> count.table
dds0 <- DESeqDataSetFromMatrix(countData=x, colData=t, design = ~ group)
print(dds0)
slotNames(dds0)
dds.norm <-  estimateSizeFactors(dds0)
sizeFactors(dds.norm)

col.sample <- c("Lgl_RNAi.24H"="orange","Lgl_RNAi.96H"="red", "WT.24H"="gray77", "WT.96H"="gray29")
t$color <- col.sample[as.vector(t$group)]
col.pheno.selected <- t$color

#Diagnostic Plots:
par(mfrow=c(2,2),cex.lab=0.7)
boxplot(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected, cex.axis=0.7,
las=1, xlab="log2(counts)", horizontal=TRUE, main="Raw counts")
boxplot(log2(counts(dds.norm, normalized=TRUE)+epsilon),  col=col.pheno.selected, cex.axis=0.7,
las=1, xlab="log2(normalized counts)", horizontal=TRUE, main="Normalized counts")
plotDensity(log2(counts(dds.norm)+epsilon),  col=col.pheno.selected,
xlab="log2(counts)", main="Density plot for Raw counts", cex.lab=0.7, panel.first=grid())
plotDensity(log2(counts(dds.norm, normalized=TRUE)+epsilon), col=col.pheno.selected,
xlab="log2(normalized counts)", main="Density plot for Normalized counts", cex.lab=0.7, panel.first=grid())

norm.counts <- counts(dds.norm, normalized=TRUE)
mean.counts <- rowMeans(norm.counts)
variance.counts <- apply(norm.counts, 1, var)
#sum(mean.counts==0) # Number of completely undetected genes
norm.counts.stats <- data.frame(min=apply(norm.counts, 2, min),
				mean=apply(norm.counts, 2, mean),
				median=apply(norm.counts, 2, median),
				max=apply(norm.counts, 2, max),
				zeros=apply(norm.counts==0, 2, sum),
				percent.zeros=100*apply(norm.counts==0, 2, sum)/nrow(norm.counts),
				perc05=apply(norm.counts, 2, quantile, 0.05),
				perc10=apply(norm.counts, 2, quantile, 0.10),
				perc90=apply(norm.counts, 2, quantile, 0.90),
				perc95=apply(norm.counts, 2, quantile, 0.95))
kable(norm.counts.stats)

par(mfrow=c(1,1))
mean.var.col <- densCols(x=log2(mean.counts), y=log2(variance.counts))
plot(x=log2(mean.counts),
			y=log2(variance.counts),
			pch=16,
			cex=0.5,
			col=mean.var.col, main="Mean-variance relationship",
			xlab="Mean log2(normalized counts) per gene",
			ylab="Variance of log2(normalized counts)",
			panel.first = grid())
			abline(a=0, b=1, col="brown")

p <- 1/6 # the probability of success
n <- 1   # target for number of successful trials
#The density function
plot(0:10, dnbinom(0:10, n, p), type="h", col="blue", lwd=4)
dnbinom(0, n , p)
sum(dnbinom(0:5, n , p)) # == pnbinom(5, 1, p)
sum(dnbinom(0:10, n , p)) # == pnbinom(10, 1, p)
1-sum(dnbinom(0:10, n , p)) # == 1 - pnbinom(10, 1, p)
n <- 2
plot(0:30, dnbinom(0:30, n, p), type="h", col="blue", lwd=2,
main="Negative binomial density",
ylab="P(x; n,p)",
xlab=paste("x = number of failures before", n, "successes"))
#Expected value
q <- 1-p
(ev <- n*q/p)
abline(v=ev, col="darkgreen", lwd=2)
#Variance
(v <- n*q/p^2)
arrows(x0=ev-sqrt(v), y0 = 0.04, x1=ev+sqrt(v), y1=0.04, col="brown",lwd=2, code=3, , length=0.2, angle=20)
n <- 10
p <- 1/6
q <- 1-p
mu <- n*q/p
all(dnbinom(0:100, mu=mu, size=n) == dnbinom(0:100, size=n, prob=p))

#Performing estimation of dispersion parameter
dds.disp <- estimateDispersions(dds.norm)

#Diagnostic plot which shows the mean of normalized counts (X-axis) and dispersion estimate for each genes
plotDispEsts(dds.disp)
alpha <- 0.0001
wald.test <- nbinomWaldTest(dds.disp)
res.DESeq2 <- results(wald.test, alpha=alpha, pAdjustMethod="BH")

#Perform the hierarchical clustering with a distance based on Pearson-correlation coefficient and average linkage clustering as agglomeration criteria
#Gene names are selected based on FDR (1%)
gene.kept <- rownames(res.DESeq2)[res.DESeq2$padj <= alpha & !is.na(res.DESeq2$padj)]
#Normalized counts for genes of interest are retrieved
count.table.kept <- log2(count.table + epsilon)[gene.kept, ]
dim(count.table.kept)
#The following heatmap is provided as Fig.1e:
heatmap.2(as.matrix(count.table.kept),
				scale="row",
				hclust=function(x) hclust(x,method="average"),
				distfun=function(x) as.dist((1-cor(t(x)))/2),
				trace="none",
				density="none",
				labRow="",
				cexCol=0.7)

dds0<-DESeqDataSetFromMatrix(countData=x, colData=t, design = ~ group)
dds <- DESeq(dds0, betaPrior=FALSE)
#The following PCA Plot is provided as Fig.1d:
plotPCA(varianceStabilizingTransformation(dds), intgroup=c("Genotype")) + theme + theme(legend.position = "top")

#EDGER FOR DIFFERENTIAL GENE EXPRESSION ANALYSIS
targets <- read.csv("~/sample_info.csv")
#Under column named PC1 in sample_info.csv, 96h-LglRNAi replicates are labeled "yes" while other samples and their replicates are labeled "no".
#This sets up the contrast observed across PC1 in Fig.1d: 96h-tjTS>lglRNAi vs all other samples (24h-tjTS>GFP + 96h-tjTS>GFP + 24h-tjTS>lglRNAi).
group <- factor(targets$PC1)
t <- cbind(targets,group=group)
x <- read.csv("~/featurecounts_bulkSeq.csv", row.names="Geneid")

genetable=data.frame(gene_id=rownames(x))
y <- DGEList(counts=x, group=group, genes=genetable)
cperm <- cpm(y)
Summary <- summary(cperm)
countcheck <- cperm > 1
keep <- which(rowSums(countcheck) >= 2)
y <- y[keep, ]
logcounts <- cpm(y,log=TRUE)

y <- calcNormFactors (y,method='TMM')
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)

y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit <- glmQLFit(y,design)
lrt<-glmLRT(fit)

my.contrasts <- makeContrasts(PC1=yes-no, levels=design)

qlf.PC1 <- glmQLFTest(fit, contrast=my.contrasts[,"PC1"])
qlf.PC1 <- glmTreat(fit, contrast=my.contrasts[,"PC1"], lfc=log2(1.2))

de1 <- decideTestsDGE(qlf.PC1, adjust.method="BH", p.value=0.05)
tab1 <- topTags(qlf.PC1, n=Inf, sort.by="PValue", p.value=0.05)

#The following file is provided as Supplementary Data 1 in the manuscript:
write.csv(tab1, file="~/DGEgenes_LglRNAi96HvsOthers(PC1).csv")
