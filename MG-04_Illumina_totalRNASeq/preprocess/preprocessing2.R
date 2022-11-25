library(knitr)
library(RColorBrewer)
library(SummarizedExperiment)
library(geneplotter)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(edgeR)
library(corpcor)  
library(M3C)
library(calibrate)
library(EnhancedVolcano) 
library(sva) 
library(pvclust) # Hierarchical clustering with p-values
library(dendextend) # fancy plotting dendograms
library(caret)
library(psych)
library(sigclust2)

setwd("/home/bscuser/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess")


#Creating SummarizedExperiment object
se <- readRDS("patients_and_controls")

# Quality assessment and normalization

## Preprocessing data
dge <- DGEList(counts = assays(se)$counts, genes = as.data.frame(mcols(se)), group = se$Case.Control)

# Computing counts per million
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM[5:10, 5:10]

dge$counts[5:10, 5:10]
dge$samples$lib.size

## Examining the sequencing depth
ord <- order(dge$samples$lib.size)
dge$samples$lib.size[ord]


# Barplot patients/controls
barplot(dge$samples$lib.size[ord]/1e+06, las = 1, main="Barplot - library size", ylab = "Million of counts", xlab = "Samples", col = c("blue","red")[(se$Case.Control[ord] == "case") +1], border = NA)
legend("topleft", c("cases","control"), fill = c("red","blue"), inset = 0.01)

# QQplot
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
qqnorm(sampledepth)
qqline(sampledepth)

# Distribution of expression levels among genes
avgexp <- rowMeans(assays(se)$logCPM)
hist(avgexp, xlab="log2 CPM", main="Distribution of the genes average expression", las=1)
abline(v=1, col="red", lwd=2)


# Filtering lowly expressed genes
nsamples <- length(se$Process.ID)
sample_cutoff <- 0.2
nsamples_cutoff <- sample_cutoff*nsamples
nsamples_cutoff

logcpm_cutoff <- 1
mask <- rowSums(assays(se)$logCPM <= logcpm_cutoff) >= nsamples_cutoff

dim(se)
se.filt <- se[!mask, ]
dim(se.filt)
dge.filt <- dge[!mask, ]
dim(dge.filt)
kept_genes2b <- dim(se.filt)[1]

# Between sample normalization
dge.filt.norm <- calcNormFactors(dge.filt)
assays(se.filt)$logCPM.norm <- cpm(dge.filt.norm, log = TRUE, prior.count = 3, normalized.lib.sizes = TRUE)

# Plots: in-sample normalized only, filtered, filtered and normalized between samples

multidensity(as.list(as.data.frame(assays(se[, se$Case.Control == "case"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM.norm)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)


# PCA based on normalized counts(normalized between samples)
pca_norm <- prcomp(as.data.frame(t(assays(se.filt)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm, var.axes=FALSE, groups=se.filt$Case.Control, labels = se.filt$Sample.ID) + ggtitle("PCA based on normalized data and colored by type")

logCPM <- cpm(dge.filt.norm, log=TRUE, prior.count=3)


## Finding highly variable genes

v <- apply(assay(se.filt), 1, var)
df <- as.data.frame(v)
colnames(df) <- 'Variance'


sel <- order(v, decreasing=TRUE)[1:1100]
se.sel <- se.filt[sel]
dge.sel <- DGEList(counts = assays(se.sel)$counts, genes = as.data.frame(mcols(se.sel)), group = se.sel$Case.Control)


#PCA based only on genes with highest variance

pca_norm_sel <- prcomp(as.data.frame(t(assays(se.sel)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_sel, var.axes=FALSE, groups=se.sel$Case.Control, labels = se.sel$Sample.ID) + ggtitle("PCA genes with higher variability")

logCPM_sel <- cpm(dge.sel, log=TRUE, prior.count=3)

counts.sel <- as.data.frame(dge.sel$counts)
row.names(counts.sel) <- rownames(dge.sel$
                                    counts)
write.csv(counts.sel, "counts_sel.csv")

##Finding correlated genes, creating correlation matrix

c <- corr.test(t(assay(se.sel)), method="spearman")
corr.df <- as.data.frame(c$r)
corr.df.p <- as.data.frame(c$p)
corr.df.t <- as.data.frame(c$t)
write.csv(corr.df, "correlation_matrix.csv")
write.csv(corr.df.p, "correlation_p.csv")
write.csv(corr.df.t, "correlation_t.csv")





