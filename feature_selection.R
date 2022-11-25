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
setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq")
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

gene.lst <- c(gene.lst,n.sel)
gene.lst <- unlist(gsub("\\.","-",gene.lst), use.names = FALSE)
gene.counts <- read.csv("../counts_sel.csv", header=TRUE, row.names=1) + 1
gene.lst.count <- gene.counts[gene.lst,]
rownames(gene.lst.count) <- gene.lst
l2 <- lapply(gene.lst.count,log2)
l2 <- as.data.frame(lapply(l2,setNames,nm=gene.lst))

#Alternative: Using selected genes from cliques (calculated in python using networkx)
gene.lst <- colnames(read.csv("correlation/cliq_gene_sel.csv"))
gene.counts <- read.csv("counts/counts_sel_overmean.csv", header=TRUE, row.names=1) + 1
gene.lst.count <- gene.counts[gene.lst,]
l2 <- lapply(gene.lst.count,log2)
l2 <- as.data.frame(lapply(l2,setNames,nm=gene.lst))

#Finding fold change between control and each patient
v_ctrl <- l2[grep("HU",colnames(l2))]
l_ctrl <- l2[grep("HD",colnames(l2))]
patients <- l2[-(grep("H",colnames(l2)))]

fold.change <- function(p,ctrl) {
  f.c <- c()
  for (i in 1:length(names(p))) {
    f.c[i] <- p[i]-ctrl[i]
  } 
  return(f.c)
}

p_mean <- apply(patients, 1, mean)
v_mean <- apply(v_ctrl, 1, mean)
l_mean <- apply(l_ctrl, 1, mean)
v_d <- as.data.frame(fold.change(p_mean, v_mean))
l_d <- as.data.frame(fold.change(p_mean, l_mean))
rownames(v_d) <- rownames(patients)
rownames(l_d) <- rownames(patients)
colnames(v_d) <- "Difference"
colnames(l_d) <- "Difference"

dif.lst.v <- list(1:36)
dif.lst.l <- list(1:36)

for (i in 1:length(colnames(v_d))) {
  dif.lst.v[i] <- list(subset(v_d[i],abs(v_d[[i]])>sort(abs(v_d[[i]]),decreasing = TRUE)[300]))
}

for (i in 1:length(colnames(l_d))) {
  dif.lst.l[i] <- list(subset(l_d[i],abs(l_d[[i]])>sort(abs(l_d[[i]]),decreasing = TRUE)[300]))
}

dif.lst.v <- as.data.frame(dif.lst.v)
dif.lst.l <- as.data.frame(dif.lst.l)

gene.names <- list(rownames(dif.lst.v),rownames(dif.lst.l))


names(gene.names) <- c("Venous", "Lymphatic")
capture.output(gene.names, file="sel_genes.txt")
saveRDS(gene.names, file="sel_genes")
