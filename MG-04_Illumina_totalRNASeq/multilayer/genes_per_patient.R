
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
library(RColorBrewer)


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
others <- c("VM119",
            "VM103",
            "VM081",
            "VM093",
            "VM090",
            "VM124",
            "VM082",
            "VM055",
            "VM068",
            "VM085",
            "VM089",
            "VM125",
            "VM092",
            "VM054"
)

other.con <- c("strange case",
               "capillary malformation", 
               "capillary malformation",
               "eccrine angiomatous hamartoma",
               "fibroadipose vascular anomaly (FAVA): lymphatic component",
               "fibroadipose vascular anomaly (FAVA): lymphatic component",
               "glomuvenous malformation",
               "infantile hemangioma",
               "infantile hemangioma",
               "venolymphatic malformation",
               "venolymphatic malformation",
               "venolymphatic malformation",
               "venous malformation or hemangioma?",
               "venous malformation: hiperplasia"
)

conditions <- c(rep("venous malformation",length(v_known)),other.con,rep("lymphatic malformation",length(l_known)))

names(conditions) <- c(v_known,others,l_known)
conditions <- conditions[order(names(conditions))]

genes_per_patient_list.l <- list()
genes_per_patient_list.v <- list()

for(k in 0:10){
  genes_per_patient.v <- list(splited_patients.v)
  for(i in 1:length(splited_patients.v)){
    genes_interesantes.v <- names(which(table(splited_patients.v[[i]])>=1))
    
    # Filter genes to those included in the multilayer network.
    genes_interesantes.v <- genes_interesantes.v[genes_interesantes.v %in% rownames(hamming_com)]
    matricin.v <- hamming_com[genes_interesantes.v,genes_interesantes.v]
    
    
    # For each value of theta between 0 and 10, we create a vector (suc), where for each gene of the multilayer we calculate the
    # number of other multilayer genes that are below a maximum k value of theta -distance in the tree-.
    names_matricin.v <- colnames(matricin.v)
    suc.v <- vector("numeric",length = length(names_matricin.v))
    names(suc.v) <- names_matricin.v
    for(j in 1:ncol(matricin.v)){
      leng_matricin.v <- length(matricin.v[,j][matricin.v[,j]<=k]) ### k (in our case, a representation of theta) is maximum distance allowed in the clustering plot. (Sauron's eye)
      suc.v[j] <- leng_matricin.v
    }
    suc.v <- suc.v[suc.v>1] ###Filetr genes that have no patient associated partners
    genes_per_patient.v[[i]] <- suc.v - 1 ### As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
  }
  names(genes_per_patient.v) <- names(splited_patients.v)
  genes_per_patient_list.v[[k+1]] <- genes_per_patient.v
  names(genes_per_patient_list.v)[[k+1]] <- as.character(k)
  message(paste0("theta= ",k))
}
  
for(k in 0:10){
  genes_per_patient.l <- list(splited_patients.l)
  for(i in 1:length(splited_patients.l)){
    genes_interesantes.l <- names(which(table(splited_patients.l[[i]])>=1))
    
    # Filter genes to those included in the multilayer network.
    genes_interesantes.l <- genes_interesantes.l[genes_interesantes.l %in% rownames(hamming_com)]
    matricin.l <- hamming_com[genes_interesantes.v,genes_interesantes.l]
    
    
    # For each value of theta between 0 and 10, we create a vector (suc), where for each gene of the multilayer we calculate the
    # number of other multilayer genes that are below a maximum k value of theta -distance in the tree-.
    names_matricin.l <- colnames(matricin.l)
    suc.l <- vector("numeric",length = length(names_matricin.l))
    names(suc.l) <- names_matricin.l
    for(j in 1:ncol(matricin.l)){
      leng_matricin.l <- length(matricin.l[,j][matricin.l[,j]<=k]) ### k (in our case, a representation of theta) is maximum distance allowed in the clustering plot. (Sauron's eye)
      suc.l[j] <- leng_matricin.l
    }
    suc.l <- suc.l[suc.l>1] ###Filter genes that have no patient associated partners
    genes_per_patient.l[[i]] <- suc.l - 1 ### As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
  }
  names(genes_per_patient.l) <- names(splited_patients.l)
  genes_per_patient_list.l[[k+1]] <- genes_per_patient.l
  names(genes_per_patient_list.l)[[k+1]] <- as.character(k)
  message(paste0("theta= ",k))
}

names(genes_per_patient.l) <- paste0(names(genes_per_patient.l), ".l")
names(genes_per_patient.v) <- paste0(names(genes_per_patient.v), ".v")

v.and.l <- c(genes_per_patient.v, genes_per_patient.l)
dif.lst.v <- list(seq(1,36))
dif.lst.l <- list(seq(1,36))
int.lst <- list(seq(1,36))
union.lst <- list(seq(1,36))

for (j in 1:(length(v.and.l)/2)) {
  dif.lst.v[[j]] <- setdiff(names(v.and.l[[j]]), names(v.and.l[[j+36]]))
  dif.lst.l[[j]] <- setdiff(names(v.and.l[[j+36]]), names(v.and.l[[j]]))
  int.lst[[j]] <- intersect(names(v.and.l[[j+36]]), names(v.and.l[[j]]))
  }



names(dif.lst.v) <- rownames(patient_matrix.l)
names(dif.lst.l) <- rownames(patient_matrix.l)
names(int.lst) <- rownames(patient_matrix.l)



pct <- c()

for (e in 1:length(names(dif.lst.v))) {
  both <- length(dif.lst.v[[e]]) + length(dif.lst.l[[e]])
  pct[e] <- length(dif.lst.v[[e]])/both-length(dif.lst.l[[e]])/both
}
names(pct) <- rownames(patient_matrix.l)



ord <- order(pct)
odd <- seq(1,36,by=2)
x <- data.frame(pct,1)
pct <- pct[order(names(pct))]
pct <- pct[order(conditions)]
conditions <- conditions[order(conditions)]

df <- data.frame(pct,conditions)
df$conditions <- as.factor(df$conditions)
df$conditions_new <- as.numeric(df$conditions)
df <- df[order(df$pct),]
height <- rep(1,length(df$pct))
df <- data.frame(df,height)

palette(brewer.pal(n = 11, name = "Paired"))
plot(df$pct, df$height, type="o",pch=19,axes=F,xlab="",ylab="", col=df$conditions_new)
text(x[ord[odd],], rownames(x[ord[odd],]),xpd=T, pos=3, srt=90,cex=1, offset=1.7)
text(x[-c(ord[odd]),], rownames(x[-c(ord[odd]),]),xpd=T, pos=1, srt=90,cex=1, offset=1.7)


