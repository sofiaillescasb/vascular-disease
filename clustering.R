library(sigclust2)
library(jaccard)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("~/Vascular_Disease")
g.lst <- readRDS("~/Vascular_Disease/sel_genes_per_patient")
patient_matrix <- t(data.frame(lapply(g.lst,function(x) as.integer(unique(unlist(g.lst)) %in% x))))
colnames(patient_matrix) <- unique(unlist(g.lst))


vfun <- function(x, y) {1 - jaccard(x, y)}
mfun <- function(x) {
  as.dist(outer(split(x, f = row(x)), split(x, f = row(x)),
                Vectorize(vfun)))
}

set.seed(2022)
res_shc <- shc(patient_matrix, matmet=mfun, linkage="ward.D2", n_sim = 1000)
res_shc$hc_dat$labels <- rownames(patient_matrix)

png(file="plots/hc/hc_pre_multi.png", width =465, height = 225, units = "mm", res=300)

plot(res_shc,alpha=0.5,ci_emp=T,use_labs = TRUE)

dev.off()

#Mapping gene IDs

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(t(patient_matrix)), keytype="SYMBOL", column="ENTREZID", multiVals="first")
head(geneSymbols)
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds]

g.lst.eid <- lapply(g.lst, function(x) found_genes[x])


# using CmmD to generate base networks from files containing data extracted from databases

redes <- paste0('gene_multilayer_network-master/networks/',list.files('gene_multilayer_network-master/networks/'))

#The next line was run in a separate bash terminal, because RStudio doesn't source from bashrc
#com_results <- CmmD::CmmD(input_layers = redes,resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = 'hamming',threads = 7,destfile_community_analysis = 'Com_Out/')

ham_dist <- read.csv("communities/hamming_distance_multilayer_network.tsv", sep ='')
sel.ham <- lapply(g.lst.eid, function(x) ham_dist[x,])
