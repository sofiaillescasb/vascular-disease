library(caret)
library(psych)
library(RCy3)
library(igraph)

setwd("~/Vascular_Disease")

correlations <- read.csv("correlation/correlation_matrix.csv", header=TRUE, row.names=1)
correlations.2 <- correlations
correlated.genes <- correlations

#Adjacency matrix
correlations.2[abs(correlations.2) > 0.9] <- 1
correlations.2[abs(correlations.2) <= 0.9] <- 0
correlations.2 <- as.matrix(correlations.2)


#Removing self correlations
diag(correlations.2) <- 0

write.csv(correlations.2, "~/Vascular_Disease/correlation/correlations_2.csv")
corr.graph <- igraph::graph_from_adjacency_matrix(as.matrix(correlations.2))
components(corr.graph)
createNetworkFromIgraph(corr.graph)
edge.list <- as_edgelist(corr.graph, names = TRUE)
write.csv(edge.list, "~/Vascular_Disease/edges.csv")

#Returns columns to reduce, used here to separate selected and discarded genes into separate objects
highCorr <- findCorrelation(correlations, cutoff = .9, names = FALSE)
adj_matr <- correlations[, -highCorr]
filtered <- correlations[, highCorr]
write.csv(adj_matr, "~/Vascular_Disease/correlation/adjacency_matrix.csv")
write.csv(filtered, "~/Vascular_Disease/correlation/reduced_genes.csv")
correlated.genes[abs(correlated.genes) < 0.9 | correlated.genes == 1] <- "NA"
write.csv(correlated.genes, "~/Vascular_Disease/correlation/correlated_genes.csv")

colnames(adj_matr)
colnames(filtered)

#List of selected genes:
#All uncorrelated genes and the gene with highest betweenness centrality of each group of correlated genes

corr.bet <- estimate_betweenness(corr.graph, directed = FALSE, cutoff = -1)
comps <- split(names(V(corr.graph)), components(corr.graph)$membership)

write.csv(corr.bet, "~/Vascular_Disease/correlation/bet.csv")

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

gene.lst <- c(gene.lst,n.sel)
gene.lst <- unlist(gsub("\\.","-",gene.lst), use.names = FALSE)
gene.counts <- read.csv("counts/counts_sel_overmean.csv", header=TRUE, row.names=1) + 1
gene.lst.count <- gene.counts[gene.lst,]
rownames(gene.lst.count) <- gene.lst
l2 <- lapply(gene.lst.count,log2)
l2 <- as.data.frame(lapply(l2,setNames,nm=gene.lst))

#Finding fold change between control and each patient
ctrl <- l2[1]
patients <- l2[-1]

fold.change <- function(p) {
  p-ctrl
}

d <- as.data.frame(lapply(patients, fold.change))
colnames(d) <- colnames(l2[-1])
dif.lst <- list(1:36)

for (i in 1:length(colnames(d))) {
  dif.lst[i] <- list(subset(d[i],abs(d[[i]])>sort(abs(d[[i]]),decreasing = TRUE)[300]))
}

names(dif.lst) <- colnames(l2[-1])
gene.names <- c()

for (j in 1:length(dif.lst)) {
  gene.names <- c(gene.names,list(rownames(dif.lst[[j]])))
  
}

names(gene.names) <- colnames(l2[-1])
capture.output(gene.names, file="~/Vascular_Disease/sel_genes_per_patient.txt")
saveRDS(gene.names, file="~/Vascular_Disease/sel_genes_per_patient")
