library(sigclust2)
library(jaccard)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gprofiler2)
library(dendroextras)
library(dendextend)
library(circlize)

setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/feature_selection")
g.lst <- readRDS("sel_genes_v")
patient_lst <- lapply(g.lst,function(x) as.integer(unique(unlist(g.lst)) %in% x))
patient_matrix <- lapply(patient_lst, function(x) as.data.frame(x))

for (i in 1:length(patient_matrix)) {
  rownames(patient_matrix[[i]]) <- unique(unlist(g.lst))
}

patient_matrix <- lapply(patient_matrix, function(x) as.data.frame(x))
patient_matrix <- as.data.frame(patient_matrix)
colnames(patient_matrix) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                              "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                              "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                              "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                              "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                              "VM125", "VM127")
patient_matrix <- t(patient_matrix)
vfun <- function(x, y) {1 - jaccard(x, y)}
mfun <- function(x) {
  as.dist(outer(split(x, f = row(x)), split(x, f = row(x)),
                Vectorize(vfun)))
}

set.seed(2022)
res_shc <- shc(patient_matrix, matmet=mfun, linkage="ward.D2", n_sim = 1000)
res_shc$hc_dat$labels <- rownames(patient_matrix)
saveRDS(res_shc, "res_shc")
png(file="hc_pre_multi.png", width =465, height = 225, units = "mm", res=300)
plot(res_shc,alpha=0.5,ci_emp=T,use_labs = TRUE)
dev.off()

#Mapping gene IDs

geneSymbols <- mapIds(org.Hs.eg.db, keys=rownames(t(patient_matrix)), keytype="SYMBOL", column="ENTREZID", multiVals="first")
head(geneSymbols)
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds]

g.lst.eid <- lapply(g.lst, function(x) found_genes[x])


# using CmmD to generate base networks from files containing data extracted from databases

#redes <- paste0('gene_multilayer_network-master/networks/',list.files('gene_multilayer_network-master/networks/'))

#The next line was run in a separate bash terminal, because RStudio doesn't source from bashrc
#com_results <- CmmD::CmmD(input_layers = redes,resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = 'hamming',threads = 7,destfile_community_analysis = 'Com_Out/')

ham_dist <- read.csv("communities/hamming_distance_multilayer_network.tsv", sep ='')

ifun <- function(g) {
  inx <- c()
  for (k in 1:(length(g))) {
      inx <- c(inx,which(rownames(ham_dist)==g[k]))
  }
  return(inx)
}

index.lst <- lapply(g.lst.eid, function(x) ifun(x))

#Making a list of the patients' dataframes
df.lst <- lapply(index.lst, function(x) ham_dist[x,x])

cfun <- function(d,n) {
  df.filt <- apply(d, 1, function(x) which(x!=24))
}

df.filt.lst <- lapply(df.lst, cfun)
df.filt.lst2 <- df.filt.lst

#Removing genes that only form communities with themselves

for (m in 1:length(df.filt.lst2)) {
  len <- sapply(df.filt.lst[[m]], length)
  df.filt.lst2[[m]] <- df.filt.lst[[m]][order(len)]
  if (length(df.filt.lst2[[m]][[1]])==1) {
    df.filt.lst2[[m]] <- df.filt.lst2[[m]][-(1)]
  }
}

post.mn.patients <- lapply(df.filt.lst2, function(x) names(x))
post.mn.patients.matrix <- t(data.frame(lapply(post.mn.patients,function(x) as.integer(unique(unlist(post.mn.patients)) %in% x))))
saveRDS(df.lst, file="~/Vascular_Disease/sel_hamming")
saveRDS(post.mn.patients, file="~/Vascular_Disease/filt_genes_per_patient")


set.seed(2022)
res_shc2 <- shc(post.mn.patients.matrix, matmet=mfun, linkage="ward.D2", n_sim = 1000)
res_shc2$hc_dat$labels <- rownames(post.mn.patients.matrix)

png(file="~/Vascular_Disease/plots/hc/hc_post_multi.png", width =465, height = 225, units = "mm", res=300)

plot(res_shc2,alpha=0.5,ci_emp=T,use_labs = TRUE)

dev.off()


