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

#DIFFERENTIAL EXPRESSION

gene.lst <- c(gene.lst,n.sel)
gene.lst <- unlist(gsub("\\.","-",gene.lst), use.names = FALSE)
se.sel2 <- se.sel[gene.lst,]
assay(se.sel2) <- as.matrix(assay(se.sel2))
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
genes <- rownames(resOrdered[1:166,])

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
genes <- rownames(resOrdered_b[1:10,])
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


