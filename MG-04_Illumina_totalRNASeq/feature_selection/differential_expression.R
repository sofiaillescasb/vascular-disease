setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/")

library(DESeq2)
library(caret)
library(psych)
library(RCy3)
library(igraph)


correlations <- read.csv("preprocess/correlation_matrix.csv", header=TRUE, row.names=1)
correlations.2 <- correlations
correlated.genes <- correlations

#Adjacency matrix
correlations.2[abs(correlations.2) > 0.9] <- 1
correlations.2[abs(correlations.2) <= 0.9] <- 0
correlations.2 <- as.matrix(correlations.2)


#Removing self correlations
diag(correlations.2) <- 0
write.csv(correlations.2, "feature_selection/correlations_2.csv")
corr.graph <- igraph::graph_from_adjacency_matrix(as.matrix(correlations.2))
components(corr.graph)
createNetworkFromIgraph(corr.graph)
edge.list <- as_edgelist(corr.graph, names = TRUE)
write.csv(edge.list, "feature_selection/edges.csv")

#Returns columns to reduce, used here to separate selected and discarded genes into separate objects
highCorr <- findCorrelation(correlations, cutoff = .9, names = FALSE)
adj_matr <- correlations[, -highCorr]
filtered <- correlations[, highCorr]
write.csv(adj_matr, "feature_selection/adjacency_matrix.csv")
write.csv(filtered, "feature_selection/reduced_genes.csv")
correlated.genes[abs(correlated.genes) < 0.9 | correlated.genes == 1] <- "NA"
write.csv(correlated.genes, "feature_selection/correlated_genes.csv")

colnames(adj_matr)
colnames(filtered)

#List of selected genes:
#All uncorrelated genes and the gene with highest betweenness centrality of each group of correlated genes

corr.bet <- estimate_betweenness(corr.graph, directed = FALSE, cutoff = -1)
comps <- split(names(V(corr.graph)), components(corr.graph)$membership)

write.csv(corr.bet, "bet.csv")

corr.comp <- components(corr.graph)
corr.tab <- table(corr.comp$membership)
corr <- corr.tab[corr.tab>1]
uncorr <- corr.tab[corr.tab<2]
corr.groups <- comps[names(corr)]

#Final gene list
gene.lst <- unname(comps[names(uncorr)])

#Selecting the correlated genes with higher betweenness centrality
n.sel<- c()

for (g in corr.groups){
  if (length(g) != 2){
    n.sel <- c(n.sel, names(which.max(corr.bet[g])))
  }
  else {
    n.sel <- c(n.sel, g[1])
  }
}

se.sel <- readRDS("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/se_sel")
se.filt <- se.sel

###################### DEAnalysis 1 ############################################
# Without correcting
se.filt$Case.Control <- relevel(se.filt$Case.Control, ref="control")
mod <- model.matrix(~ se.filt$Case.Control, colData(se.filt))
mod0 <- model.matrix(~ 1, colData(se.filt))
pv <- f.pvalue(assays(se.filt)$logCPM.norm, mod, mod0)
sum(p.adjust(pv, method="fdr") < 0.05)
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



#DIFFERENTIAL EXPRESSION 2

gene.lst <- c(gene.lst,n.sel)
gene.lst <- unlist(gsub("\\.","-",gene.lst), use.names = FALSE)
se.sel2 <- se.sel
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

resLFC <- lfcShrink(dds, coef="Case.Control_control_vs_case", type="apeglm")
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
genes <- rownames(resOrdered_b[1:66,])
genes

# res05_l <- results(dds, contrast=c("group","case", "l_control"), alpha=0.05)
# summary(res05_l)
# sum(res05_l$padj < 0.05, na.rm=TRUE)
# res05_l$padj[res05_l$padj < 0.05]
# 
# resLFC_v <- lfcShrink(dds, coef="group_v_control_vs_case", type="apeglm")
# resLFC_v
# res_vOrdered <- res05_v[order(res05_v$padj),]
# genes_v <- rownames(res_vOrdered[1:149,])


