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

setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/")


#Importing data as summarizedexperiment object
se <- readRDS("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/patients_and_controls")

# Data exploration
metadata <- as.data.frame(colData(se))
fplot <- ggplot(as.data.frame(metadata))
head(metadata)

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
#barplot(dge$samples$lib.size[ord]/1e+06, las = 1, main="Barplot - library size", ylab = "Million of counts", xlab = "Samples", col = c("blue","red")[(se$Case.Control[ord] == "case") +1], border = NA)
#legend("topleft", c("cases","control"), fill = c("red","blue"), inset = 0.01)

# QQplot
#sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
#qqnorm(sampledepth)
#qqline(sampledepth)

# Distribution of expression levels among genes
#avgexp <- rowMeans(assays(se)$logCPM)
#hist(avgexp, xlab="log2 CPM", main="Distribution of the genes average expression", las=1)
#abline(v=1, col="red", lwd=2)


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

saveRDS(se.filt, "~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/se_filt")

# Plots: in-sample normalized only, filtered, filtered and normalized between samples

#multidensity(as.list(as.data.frame(assays(se[, se$Case.Control == "case"])$logCPM)),
             #xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

#multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM)),
 #            xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

#multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM.norm)),
 #            xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)


# PCA based on normalized counts(normalized between samples)
pca_norm <- prcomp(as.data.frame(t(assays(se.filt)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm, var.axes=FALSE, groups=se.filt$Case.Control, labels = se.filt$Sample.ID) + ggtitle("PCA based on normalized data and colored by type")

###################### DEAnalysis 1 ############################################
# Without correcting for batch

se.filt$Case.Control <- relevel(factor(se.filt$Case.Control), ref="control")
mod <- model.matrix(~ se.filt$Case.Control, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))
pv <- f.pvalue(assays(se.filt)$logCPM.norm, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)
names(pv[p.adjust(pv, method="fdr") < 0.05])
# 1811 DE genes

# Distribution of expected p-values
par(mfrow=c(1, 1))
hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000))

# Correcting for the series
se.filt$Batch <- c(1:43)
se.filt$Batch[1:37] <- 1 
se.filt$Batch[38:43] <- 2 
if (length(unique(se.filt$Batch)) > 1){
  mod_b <- model.matrix(~Case.Control + Batch, colData(se.filt))
  mod0_b <- model.matrix(~Batch, colData(se.filt))
  colnames(mod_b)
  pValues_b <- f.pvalue(assays(se.filt)$logCPM,mod_b,mod0_b)
  sum(p.adjust(pValues_b, method="fdr") < 0.05)
} else {
  print("There is only one series for this disease")
}

# Distribution of expected p-values after correcting for the series
if (length(unique(se.filt$Batch)) > 1){
  hist(pValues_b, main="Distribution of expected p-values after correcting for series", las=1)
}

# Mean-variance relationship
FDRcutoff <- 0.05
v <- voom(dge.filt.norm, mod, plot = TRUE) # Voom is applied to the normalized and filtered DGEList object
fit <- lmFit(v, mod)
fit<- eBayes(fit)
res <- decideTests(fit, p.value = FDRcutoff)
summary(res)

# Add gene metadata (gene symbol) and fetch the table of results. Print first 15 sDEgenes
rowRanges(se.filt)
#genesmd <- data.frame(chr = as.character(seqnames(rowRanges(se.filt))), symbol = rowData(se.filt)[, 1], stringsAsFactors = FALSE)
genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt, 15)
genes_m1 <- names(pv[p.adjust(pv, method="fdr") < 0.05])[order(names(pv[p.adjust(pv, method="fdr") < 0.05]))]


#DIFFERENTIAL EXPRESSION 2

se.sel2 <- se.filt
assay(se.sel2) <- as.matrix(assay(se.sel2))
se.sel2$Batch <- c(1:43)
se.sel2$Batch[1:37] <- 1 
se.sel2$Batch[38:43] <- 2 
se.sel2$Batch <- factor(se.sel2$Batch)

#Not correcting for batch
dds <- DESeqDataSet(se.sel2, design = ~ Case.Control)
dds <- DESeq(dds)
resultsNames(dds)

res05 <- results(dds, contrast=c("Case.Control","case", "control"), alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
res05$padj[res05$padj < 0.05]

resLFC <- lfcShrink(dds, coef="Case.Control_case_vs_control", type="apeglm")
resLFC
resOrdered <- res05[order(res05$pvalue),]
genes <- rownames(resOrdered[1:2287,])

#correcting for batch
dds_b <- DESeqDataSet(se.sel2, design = ~ Case.Control + Batch)
dds_b <- DESeq(dds_b)
dds_b

res05_b <- results(dds_b, name="Batch_2_vs_1", alpha=0.05)
summary(res05_b)
sum(res05_b$padj < 0.05, na.rm=TRUE)
res05_b$padj[res05_b$padj < 0.05]

resLFC_b <- lfcShrink(dds_b, coef="Batch_2_vs_1", type="apeglm")
resLFC_b
resOrdered_b <- res05_b[order(res05_b$pvalue),]
genes_b <- rownames(resOrdered_b[1:66,])
genes_b


comprson <- genes_m1[order(genes_m1)][genes_m1 %in% genes[order(genes)]]
length(comprson)
#1739 DE genes in common






