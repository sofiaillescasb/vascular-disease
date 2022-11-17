
jaccard_ind <- function(x){
  res <- matrix(data=NA,nrow=ncol(x),ncol=ncol(x))
  rownames(res) <- colnames(x)
  colnames(res) <- colnames(x)
  for(i in 1:ncol(x)){
    uno <- rownames(x)[which(x[,i]==1)]
    for(j in 1:ncol(x)){
      if(is.na(res[i,j])){
        dos <- rownames(x)[which(x[,j]==1)]
        numerador <- length(intersect(uno,dos))
        denom <- length(union(uno,dos))
        out <- numerador/denom
        res[i,j] <- out
        res[j,i] <- res[i,j]
      }
    }
  }
  attr(res,"method") <- "jaccard_ind" 
  return(as.dist(1-res)) ###Also: We have to be calculating a dissimilarity matrix with the function, so 1-x.
}


setwd('/home/bscuser/Vascular_Disease/MG-04_Illumina_totalRNASeq/multilayer')
# Load dependencies
library(sigclust2)
library(pvclust)
library(fpc)
library(jaccard)
library(CmmD)
library(readr)
library(knitr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggvenn)



sel.ham.v <- readRDS("../feature_selection/sel_hamming_v") #table generated in clustering.R 
sel.ham.l <- readRDS("../feature_selection/sel_hamming_l") #table generated in clustering.R 
hamming_com <- read.table("~/Vascular_Disease/communities/hamming_distance_multilayer_network.tsv",sep='\t') #hamming distance matrix obtained from repository, script for generation at ~/Vascular_Disease/generate_multilayer.R
colnames(hamming_com) <- rownames(hamming_com)

# Load genes associated to each patient from factor selection
tata.v <- lapply(sel.ham.v, function(x) rownames(x))
splited_patients.v <- tata.v
names(splited_patients.v) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                               "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                               "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                               "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                               "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                               "VM125", "VM127")

all_genes_possible.v <- unique(unlist(splited_patients.v,use.names=F))
all_genes_symbol.v <- mapIds(x = org.Hs.eg.db,keys = c(all_genes_possible.v), keytype = 'ENTREZID',column = 'SYMBOL')

# Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
n_genes_p_patients.v <- matrix(data= 0, nrow= length(splited_patients.v),ncol= length(all_genes_possible.v))
colnames(n_genes_p_patients.v) <- all_genes_possible.v
rownames(n_genes_p_patients.v) <- names(splited_patients.v)

for(rowi in 1:nrow(n_genes_p_patients.v)){
  n_genes_p_patients.v[rowi,] <- as.integer(colnames(n_genes_p_patients.v) %in% splited_patients.v[[rowi]])
}

patient_matrix.v <- n_genes_p_patients.v
rownames(patient_matrix.v) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                              "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                              "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                              "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                              "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                              "VM125", "VM127")

# Load genes associated to each patient from factor selection
tata.l <- lapply(sel.ham.l, function(x) rownames(x))
splited_patients.l <- tata.l
names(splited_patients.l) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                               "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                               "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                               "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                               "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                               "VM125", "VM127")

all_genes_possible.l <- unique(unlist(splited_patients.l,use.names=F))
all_genes_symbol.l <- mapIds(x = org.Hs.eg.db,keys = c(all_genes_possible.l), keytype = 'ENTREZID',column = 'SYMBOL')

# Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
n_genes_p_patients.l <- matrix(data= 0, nrow= length(splited_patients.l),ncol= length(all_genes_possible.l))
colnames(n_genes_p_patients.l) <- all_genes_possible.l
rownames(n_genes_p_patients.l) <- names(splited_patients.l)

for(rowi in 1:nrow(n_genes_p_patients.l)){
  n_genes_p_patients.l[rowi,] <- as.integer(colnames(n_genes_p_patients.l) %in% splited_patients.l[[rowi]])
}

patient_matrix.l <- n_genes_p_patients.l
rownames(patient_matrix.l) <- c("VM015", "VM024", "VM038", "VM040", "VM042", "VM043",   
                              "VM048", "VM053", "VM054", "VM055", "VM056", "VM060", "VM064",  
                              "VM066", "VM068", "VM071", "VM072", "VM073", "VM081", "VM082",   
                              "VM083", "VM085", "VM089", "VM090", "VM092", "VM093", "VM099",   
                              "VM103", "VM108", "VM110", "VM111", "VM113", "VM119", "VM124",   
                              "VM125", "VM127")

v_known <- c("VM015", "VM040", "VM042", "VM043", "VM048", "VM053", "VM056", "VM060",
             "VM064", "VM066", "VM073", "VM083", "VM099", "VM108", "VM110", "VM111",
             "VM113", "VM127")
l_known <- c("VM024", "VM038", "VM071", "VM072")

others <- rownames(patient_matrix.l)[!(rownames(patient_matrix.l)%in%v_known)]
others <- others[!(others%in%l_known)]

ven_lst <- lapply(splited_patients.v[v_known],function(x) as.integer(unique(unlist(splited_patients.v[v_known],use.names=F)) %in% x))
ven_matrix <- as.data.frame(ven_lst)
rownames(ven_matrix) <- unique(unlist(splited_patients.v[v_known],use.names=F))

lymph_lst <- lapply(splited_patients.l[l_known],function(x) as.integer(unique(unlist(splited_patients.l[l_known],use.names=F)) %in% x))
lymph_matrix <- as.data.frame(lymph_lst)
rownames(lymph_matrix) <- unique(unlist(splited_patients.l[l_known],use.names=F))

others_lst.l <- lapply(splited_patients.l[others],function(x) as.integer(unique(unlist(splited_patients.l[l_known],use.names=F)) %in% x))
others_matrix.l <- as.data.frame(others_lst.l)
rownames(others_matrix.l) <- unique(unlist(splited_patients.l[l_known],use.names=F))
others_matrix.l <- others_matrix.l*(-1)

others_lst.v <- lapply(splited_patients.v[others],function(x) as.integer(unique(unlist(splited_patients.v[v_known],use.names=F)) %in% x))
others_matrix.v <- as.data.frame(others_lst.v)
rownames(others_matrix.v) <- unique(unlist(splited_patients.v[v_known],use.names=F))

sum(others_matrix.v$VM055) + sum(others_matrix.l$VM055)
sum(others_matrix.v$VM054) + sum(others_matrix.l$VM054)

tot_v <- apply(others_matrix.v, 2, sum)
tot_l <- apply(others_matrix.l, 2, sum)

tot <- tot_v + tot_l
tot <- tot[order(tot)]
cols <- c("#0000CC", "#0000CC", "#1034A6", "#1034A6", "#412F88", "#1034A6", "#722B6A", "#722B6A", "#A2264B", "#A2264B", "#D3212D", "#D3212D", "#F62D2D", "#F62D2D")
plot(c(tot)[-c(4)], rep(1,length(tot[-c(4)])),axes=F,xlab="",ylab="",type="o",pch=19, col=cols)
text(c(tot)[-c(4)], rep(1,length(tot[-c(4)])),names(tot)[-c(4)],xpd=T, pos=1,xpd=T, srt=45,cex=0.7, offset=1)
text(c(tot)[4], rep(1,length(tot[4])),names(tot)[4],xpd=T, pos=3,xpd=T, srt=45,cex=0.7, offset=1)



