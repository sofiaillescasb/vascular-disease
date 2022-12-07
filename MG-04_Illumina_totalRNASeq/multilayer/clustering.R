library(jaccard)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(RColorBrewer)



setwd("/Users/sofiaillescas/Desktop/Vascular_Disease/vascular-disease/MG-04_Illumina_totalRNASeq/feature_selection")
g.lst.v <- readRDS("sel_genes_v")
g.lst.l <- readRDS("sel_genes_l")

meta <- read.csv("../metadata.csv",row.names = 1)

patient_lst.v <- lapply(g.lst.v,function(x) as.integer(unique(unlist(g.lst.v)) %in% x))
patient_matrix.v <- lapply(patient_lst.v, function(x) as.data.frame(x))

for (i in 1:length(patient_matrix.v)) {
  rownames(patient_matrix.v[[i]]) <- unique(unlist(g.lst.v))
}

patient_matrix.v <- lapply(patient_matrix.v, function(x) as.data.frame(x))
patient_matrix.v <- as.data.frame(patient_matrix.v)
colnames(patient_matrix.v) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                              "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                              "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                              "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                              "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                              "VM125", "VM127")
patient_matrix.v <- t(patient_matrix.v)

#Mapping gene IDs

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(t(patient_matrix.v)), keytype="SYMBOL", column="ENTREZID", multiVals="first")
head(geneSymbols)
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds]

g.lst.eid.v <- lapply(g.lst.v, function(x) found_genes[x])
saveRDS(g.lst.eid.v, file="../multilayer/genes_per_patientv2")

#Same process but using the genes that were selected as different from the lymphatic controls

patient_lst.l <- lapply(g.lst.l,function(x) as.integer(unique(unlist(g.lst.l)) %in% x))
patient_matrix.l <- lapply(patient_lst.l, function(x) as.data.frame(x))

for (i in 1:length(patient_matrix.l)) {
  rownames(patient_matrix.l[[i]]) <- unique(unlist(g.lst.l))
}

patient_matrix.l <- lapply(patient_matrix.l, function(x) as.data.frame(x))
patient_matrix.l <- as.data.frame(patient_matrix.l)
colnames(patient_matrix.l) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                              "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                              "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                              "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                              "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                              "VM125", "VM127")
patient_matrix.l <- t(patient_matrix.l)


#Mapping gene IDs

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(t(patient_matrix.l)), keytype="SYMBOL", column="ENTREZID", multiVals="first")
head(geneSymbols)
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds]

g.lst.eid.l <- lapply(g.lst.l, function(x) found_genes[x])
saveRDS(g.lst.eid.l, file="../multilayer/genes_per_patientl2")

union.lst <- list(seq(1,36))

for (j in 1:36) {
  union.lst[[j]] <- union(names(g.lst.eid.l[[j]]), names(g.lst.eid.v[[j]]))
}

names(union.lst) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                      "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                      "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                      "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                      "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                      "VM125", "VM127")

pct.v <- c()
pct.l <- c()

for (e in 1:length(names(union.lst))) {
  pct.v[e] <- length(names(g.lst.eid.v[[e]]))/length(union.lst[[e]])
  pct.l[e] <- length(names(g.lst.eid.l[[e]]))/length(union.lst[[e]])
}

names(pct.v) <- names(union.lst)
names(pct.l) <- names(union.lst)

df <- meta
df$pct.v <- pct.v
df$pct.l <- pct.l

df$conditions <- as.factor(df$conditions)
df$conditions_new <- as.numeric(df$conditions)

palette(brewer.pal(n = 11, name = "Paired"))
plot(df$pct.l, df$pct.v, type="p",pch=19,axes=T,xlab="percentage lymphatic",ylab="percentage venous", col=df$conditions_new)
text(df$pct.l, df$pct.v, rownames(df),xpd=T, pos=2,cex=0.7, offset=0.5)
abline(h=0.5)
abline(v=0.7)
legend("bottomleft",legend=unique(df$conditions),cex=0.7,fill=unique(df$conditions_new),bty = "n",y.intersp=0.6)

df$mutations <- as.factor(df$mutations)
df$mutations_new <- as.numeric(df$mutations)

palette(brewer.pal(n = 4, name = "Set2"))
plot(df$pct.l, df$pct.v, type="p",pch=19,axes=T,xlab="percentage lymphatic",ylab="percentage venous", col=df$mutations_new)
text(df$pct.l, df$pct.v, rownames(df),xpd=T, pos=2,cex=0.7, offset=0.5)
abline(h=0.5)
abline(v=0.7)
legend("bottomleft",legend=unique(df$mutations),cex=0.7,fill=unique(df$mutations_new),bty = "n",y.intersp=0.6)

df$variant <- as.factor(df$variants)
df$variant_new <- as.numeric(df$variant)

palette(brewer.pal(n = 11, name = "Paired"))
plot(df$pct.l, df$pct.v, type="p",pch=19,axes=T,xlab="percentage lymphatic",ylab="percentage venous", col=df$variant_new)
text(df$pct.l, df$pct.v, rownames(df),xpd=T, pos=2,cex=0.7, offset=0.5)
abline(h=0.5)
abline(v=0.7)
legend("bottomleft",legend=unique(df$variant),cex=0.7,fill=unique(df$variant_new),bty = "n",y.intersp=0.6)



#######################################multi###################################

ham_dist <- read.csv("../multilayer/hamming_2022.csv",row.names = 1)
ham_dist 
colnames(ham_dist) <- rownames(ham_dist)

v.genes <- unique(unlist(g.lst.eid.v))
l.genes <- unique(unlist(g.lst.eid.l))


sel.ham.v <- ham_dist[v.genes,]
sel.ham.l <- ham_dist[l.genes,]

saveRDS(sel.ham.v,file="sel_ham_v")
saveRDS(sel.ham.l,file="sel_ham_l")


