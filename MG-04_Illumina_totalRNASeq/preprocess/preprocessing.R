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
library(dendextend) # fancy plotting dendrograms
library(caret)
library(psych)
library(sigclust2)
library(ggvenn)
library(DESeq2)

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

saveRDS(se.filt, "~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/se_filt")

# Plots: in-sample normalized only, filtered, filtered and normalized between samples

multidensity(as.list(as.data.frame(assays(se[, se$Case.Control == "case"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM)),
            xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM.norm)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)


# PCA based on normalized counts(normalized between samples)
pca_norm_pre <- prcomp(as.data.frame(t(assays(se.filt)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_pre, var.axes=FALSE, groups=se.filt$Case.Control, labels = colnames(se.filt)) + ggtitle("PCA based on normalized data and colored by type")

###################### DEAnalysis ############################################

##All controls vs all patients

# Without correcting for batch
se.filt$Case.Control <- relevel(factor(se.filt$Case.Control), ref="control")
mod <- model.matrix(~ se.filt$Case.Control, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))
pv <- f.pvalue(assays(se.filt)$logCPM.norm, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)
names(pv[p.adjust(pv, method="fdr") < 0.05])
# 1811 DE genes, before mean-variance relationship

# Distribution of expected p-values
par(mfrow=c(1, 1))
hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000))

# Correcting for batch
se.filt$Batch <- c(1:43)
se.filt$Batch[1:37] <- "1" 
se.filt$Batch[38:43] <- "2" 
if (length(unique(se.filt$Batch)) > 1){
  mod_b <- model.matrix(~Case.Control + Batch, colData(se.filt))
  mod0_b <- model.matrix(~Batch, colData(se.filt))
  colnames(mod_b)
  pValues_b <- f.pvalue(assays(se.filt)$logCPM,mod_b,mod0_b)
  sum(p.adjust(pValues_b, method="fdr") < 0.05)
} else {
  print("There is only one series for this disease")
}
#No DE genes when accounting for batch number

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
# 1935 DE genes

# Add gene metadata (gene symbol) and fetch the table of results. Print first 15 sDEgenes
rowRanges(se.filt)
genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt, 15)
tt_sign <- tt[tt$adj.P.Val <= 0.05,]
head(tt_sign); nrow(tt_sign)
sum(tt$adj.P.Val < 0.05)
genes <- (tt$symbol[(tt$adj.P.Val < 0.05)])
genes_or <- genes[order(genes)]

# Q-Q plot
{qqt(fit$t[, 2], df = fit$df.prior + fit$df.residual, main = "", pch = ".", cex = 3, las = 1)
  qqline(fit$t[, 2], col = "red")}

# Volcano plot
par(xpd=T, mar=par()$mar+c(0,0,0,5))
with(tt, plot(logFC, -log10(adj.P.Val), pch = 20, main = "", xlim=c(-6,7), ylim=c(0,15))) 
with(subset(tt, logFC > 1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="green"))
with(subset(tt, logFC < -1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="orange"))
with(subset(tt, -log10(adj.P.Val) < 1.3), points(logFC, -log10(adj.P.Val), pch = 20, col="grey"))

legend("bottomright", cex = .75, inset = c(-0.32,0.75), xjust = 2, yjust = 1,pch = 20,c("sDE & overexpressed", "sDE & underexpressed", "sDE", "non-sDE"), col = c("green", "orange", "grey", "black"),bg="white", bty= "n")
with(subset(tt, -log10(adj.P.Val)>1.3 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs = symbol, cex=.7, col= "#35445b"))
par(xpd=F)
abline(h= 1.3, col = "blue", lty= 2, lwd = 1)
abline(v= c(-1,1), col = "blue", lty= 2, lwd = 1)

pca_norm_post <- prcomp(as.data.frame(t(assays(se.filt)$logCPM.norm[genes,])) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_post, var.axes=FALSE, groups=se.filt$Case.Control, labels = colnames(se.filt)) + ggtitle("PCA based on normalized data and colored by type")

## Venous controls vs venous patients
patients_v <- se.filt[,se.filt$Summary.clinic=="venous malformation"]
controls_v <- se.filt[,colnames(assay(se.filt))%in%c("HUVECS","HUVECs_2","HUVECs_3","HUVECs_4")]
se.v <- cbind(patients_v, controls_v)

se.v$Case.Control <- relevel(factor(se.v$Case.Control), ref="control")
mod <- model.matrix(~ se.v$Case.Control, colData(se.v))
mod0 <- model.matrix(~ 1, colData(se.v))
pv <- f.pvalue(assays(se.v)$logCPM.norm, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)
names(pv[p.adjust(pv, method="fdr") < 0.05])

par(mfrow=c(1, 1))
hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000))

FDRcutoff <- 0.05
v <- voom(dge.filt.norm[,colnames(se.v)], mod, plot = TRUE) # Voom is applied to the normalized and filtered DGEList object
fit <- lmFit(v, mod)
fit<- eBayes(fit)
res <- decideTests(fit, p.value = FDRcutoff)
summary(res)

genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt, 15)
tt_sign <- tt[tt$adj.P.Val <= 0.05,]
head(tt_sign); nrow(tt_sign)
sum(tt$adj.P.Val < 0.05)
genes_v <- (tt$symbol[(tt$adj.P.Val < 0.05)])
genes_v_or <- genes[order(genes)]

#Venous PCA
pca_norm_v <- prcomp(as.data.frame(t(assays(se.v)$logCPM.norm[genes_v,])) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_v, var.axes=FALSE, groups=se.v$Case.Control, labels = colnames(se.v)) + ggtitle("PCA based on venous normalized data and colored by type")


# Venous Q-Q plot
{qqt(fit$t[, 2], df = fit$df.prior + fit$df.residual, main = "", pch = ".", cex = 3, las = 1)
  qqline(fit$t[, 2], col = "red")}

# Venous volcano plot
par(xpd=T, mar=par()$mar+c(0,0,0,5))
with(tt, plot(logFC, -log10(adj.P.Val), pch = 20, main = "", xlim=c(-6,7), ylim=c(0,15))) 


with(subset(tt, logFC > 1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="green"))
with(subset(tt, logFC < -1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="orange"))
with(subset(tt, -log10(adj.P.Val) < 1.3), points(logFC, -log10(adj.P.Val), pch = 20, col="grey"))

legend("bottomright", cex = .75, inset = c(-0.32,0.75), xjust = 2, yjust = 1,pch = 20,c("sDE & overexpressed", "sDE & underexpressed", "sDE", "non-sDE"), col = c("green", "orange", "grey", "black"),bg="white", bty= "n")
with(subset(tt, -log10(adj.P.Val)>1.3 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs = symbol, cex=.7, col= "#35445b"))
par(xpd=F)
abline(h= 1.3, col = "blue", lty= 2, lwd = 1)
abline(v= c(-1,1), col = "blue", lty= 2, lwd = 1)

## Lymphatic controls vs lymphatic patients
patients_l <- se.filt[,se.filt$Summary.clinic=="lymphatic malformation"]
controls_l <- se.filt[,colnames(assay(se.filt))%in%c("HDLECs_1", "HDLECs_2", "HDLECs_3")]
se.l <- cbind(patients_l, controls_l)

se.l$Case.Control <- relevel(factor(se.l$Case.Control), ref="control")
mod <- model.matrix(~ se.l$Case.Control, colData(se.l))
mod0 <- model.matrix(~ 1, colData(se.l))
pv <- f.pvalue(assays(se.l)$logCPM.norm, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)
names(pv[p.adjust(pv, method="fdr") < 0.05])

par(mfrow=c(1, 1))
hist(pv, main="Distribution of expected p-values", las=1, ylim = c(0, 8000))

FDRcutoff <- 0.05
v <- voom(dge.filt.norm[,colnames(se.l)], mod, plot = TRUE) # Voom is applied to the normalized and filtered DGEList object
fit <- lmFit(v, mod)
fit<- eBayes(fit)
res <- decideTests(fit, p.value = FDRcutoff)
summary(res)

genesmd <- data.frame(symbol = rownames(res), stringsAsFactors = FALSE)
fit$genes <- genesmd
tt <- topTable(fit, coef = 2, n = Inf)
head(tt, 15)
tt_sign <- tt[tt$adj.P.Val <= 0.05,]
head(tt_sign); nrow(tt_sign)
sum(tt$adj.P.Val < 0.05)
genes_l <- (tt$symbol[(tt$adj.P.Val < 0.05)])
genes_l_or <- genes[order(genes)]

#Lymphatic PCA
pca_norm_l <- prcomp(as.data.frame(t(assays(se.l)$logCPM.norm[genes_v,])) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_l, var.axes=FALSE, groups=se.l$Case.Control, labels = colnames(se.l)) + ggtitle("PCA based on venous normalized data and colored by type")

# Lymphatic Q-Q plot
{qqt(fit$t[, 2], df = fit$df.prior + fit$df.residual, main = "", pch = ".", cex = 3, las = 1)
  qqline(fit$t[, 2], col = "red")}

# Lymphatic volcano plot
par(xpd=T, mar=par()$mar+c(0,0,0,5))
with(tt, plot(logFC, -log10(adj.P.Val), pch = 20, main = "", xlim=c(-6,7), ylim=c(0,15))) 


with(subset(tt, logFC > 1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="green"))
with(subset(tt, logFC < -1.0), points(logFC, -log10(adj.P.Val), pch = 20, col="orange"))
with(subset(tt, -log10(adj.P.Val) < 1.3), points(logFC, -log10(adj.P.Val), pch = 20, col="grey"))

legend("bottomright", cex = .75, inset = c(-0.32,0.75), xjust = 2, yjust = 1,pch = 20,c("sDE & overexpressed", "sDE & underexpressed", "sDE", "non-sDE"), col = c("green", "orange", "grey", "black"),bg="white", bty= "n")
with(subset(tt, -log10(adj.P.Val)>1.3 & abs(logFC)>1), textxy(logFC, -log10(adj.P.Val), labs = symbol, cex=.7, col= "#35445b"))
par(xpd=F)
abline(h= 1.3, col = "blue", lty= 2, lwd = 1)
abline(v= c(-1,1), col = "blue", lty= 2, lwd = 1)

#Compare DE genes from the three groups
genes_lst <- list(genes_l, genes_v, genes)
names(genes_lst) <- c("lymphatic", "venous", "all patients vs all controls")
saveRDS(genes_lst, file="~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/DEgenes")


ggvenn(
  genes_lst, columns=names(genes_lst),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)
