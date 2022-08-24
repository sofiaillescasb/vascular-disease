library(sigclust2)
library(jaccard)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gprofiler2)
library(dendroextras)
library(dendextend)
library(circlize)

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

png(file="plots/hc/hc_post_multi.png", width =465, height = 225, units = "mm", res=300)

plot(res_shc2,alpha=0.5,ci_emp=T,use_labs = TRUE)

dev.off()

#Dendrogram
all_genes <- unique(unlist(g.lst.eid))
all_genes2 <- c() #this vector will contain numbers to index the dataframe's columns
for (k in 1:(length(all_genes))) {
  all_genes2 <- c(all_genes2,which(rownames(ham_dist)==all_genes[k]))
}

sel <- ham_dist[all_genes2,all_genes2]
hc <- hclust(as.dist(sel,"ward.D2"))

set.seed(2022) # Set the seeds ,to generate the same dendrogram

# Generate the dendrogram object
dend3 <- as.dendrogram(hc)
 
