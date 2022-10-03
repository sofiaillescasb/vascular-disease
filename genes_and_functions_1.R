
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


setwd('Vascular_Disease/Medulloblastoma-main/')
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


sel.ham <- readRDS("~/Vascular_Disease/sel_hamming") #table generated in clustering.R 
hamming_com <- read.table("../communities/hamming_distance_multilayer_network.tsv",sep='\t') #hamming distance matrix obtained from repository, script for generation at ~/Vascular_Disease/generate_multilayer.R
colnames(hamming_com) <- rownames(hamming_com)

# Load genes associated to each patient from factor selection
tata <- lapply(sel.ham, function(x) rownames(x))
splited_patients <- tata

all_genes_possible <- unique(unlist(splited_patients,use.names=F))
all_genes_symbol <- mapIds(x = org.Hs.eg.db,keys = c(all_genes_possible), keytype = 'ENTREZID',column = 'SYMBOL')

# Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
n_genes_p_patients <- matrix(data= 0, nrow= length(splited_patients),ncol= length(all_genes_possible))
colnames(n_genes_p_patients) <- all_genes_possible
rownames(n_genes_p_patients) <- names(splited_patients)

for(rowi in 1:nrow(n_genes_p_patients)){
  n_genes_p_patients[rowi,] <- as.integer(colnames(n_genes_p_patients) %in% splited_patients[[rowi]])
}

patient_matrix <- n_genes_p_patients


vfun <- function(x, y) {1 - jaccard(x, y)}
mfun <- function(x) {
  as.dist(outer(split(x, f = row(x)), split(x, f = row(x)),
                Vectorize(vfun)))
}

set.seed(2020)
res_shc <- shc(patient_matrix, matmet=mfun, linkage="ward.D2", n_sim = 1000,alpha = .08)
res_shc$hc_dat$labels <- rownames(patient_matrix)
plot(res_shc,alpha=0.8,ci_emp=T,use_labs = TRUE)

referencia <- sigclust2::shcutree(res_shc,alpha = .08)
names(referencia) <- res_shc$hc_dat$labels

ground_truth <- data.frame(referencia)
ground_truth <- cbind(ground_truth,ground_truth[,1])
ground_truth[,1] <- rownames(ground_truth)


genes_per_patient_list <- list()
for(k in 0:10){
  genes_per_patient <- list()
  for(i in 1:length(splited_patients)){
    genes_interesantes <- names(which(table(splited_patients[[i]])>=1))
    
    # Filter genes to those included in the multilayer network.
    genes_interesantes <- genes_interesantes[genes_interesantes %in% rownames(hamming_com)]
    matricin <- hamming_com[genes_interesantes,genes_interesantes]
    
    
    # For each value of theta between 0 and 10, we create a vector (suc), where for each gene of the multilayer we calculate the
    # number of other multilayer genes that are below a maximum k value of theta -distance in the tree-.
    names_matricin <- colnames(matricin)
    suc <- vector("numeric",length = length(names_matricin))
    names(suc) <- names_matricin
    for(j in 1:ncol(matricin)){
      leng_matricin <- length(matricin[,j][matricin[,j]<=k]) ### k (in our case, a representation of theta) is maximum distance allowed in the clustering plot. (Sauron's eye)
      suc[j] <- leng_matricin
    }
    suc <- suc[suc>1] ###Filetr genes that have no patient associated partners
    genes_per_patient[[i]] <- suc - 1 ### As a value of 2 mean that only there is one more gene with the gene being analyzed, we substract 1.
  }
  names(genes_per_patient) <- names(splited_patients)
  genes_per_patient_list[[k+1]] <- genes_per_patient
  names(genes_per_patient_list)[[k+1]] <- as.character(k)
  message(paste0("theta= ",k))
}

message("Tetha based filtering finished. Calculating clustering accuracies for tetha 0 to 10 and lambda 1 to 20")

# Start a 11 x 30 matrix to be filled with the hierarchical clustering accuracy values.
final_accuracy_matrix <- matrix(0, ncol= 30, nrow= 11)
final_matthews_matrix <- matrix(0, ncol= 30, nrow= 11)
final_kk_used <- matrix(0, ncol= 30, nrow= 11)
mean_per_pair <- matrix(0, ncol= 30, nrow= 11)
rownames(final_accuracy_matrix) <- as.character(0:10)
rownames(final_kk_used) <- as.character(0:10)

u <- 5 
val <- 2  #lambda
preserve_genes_per_patient <- genes_per_patient_list[[u]] # u==k == tetha + 1
genes_per_patient <- preserve_genes_per_patient
genes_per_patient <- lapply(genes_per_patient,function(x) x[x<=val]) #val = lambda. We filter the genes that are over the lambda value tested

genes_per_patient_names <- lapply(genes_per_patient,function(x) unique(names(x)))
all_genes_possible <- unique(unlist(genes_per_patient_names,use.names=F))

# Generate a 0-1 patient x genes matrix that acts as input for pamk, jaccard_ind and hclust.
n_genes_p_patients <- matrix(data= 0, nrow= nrow(ground_truth),ncol= length(all_genes_possible))
colnames(n_genes_p_patients) <- all_genes_possible
rownames(n_genes_p_patients) <- rownames(ground_truth)

for(rowi in 1:nrow(n_genes_p_patients)){
  n_genes_p_patients[rowi,] <- as.integer(colnames(n_genes_p_patients) %in% genes_per_patient_names[[rowi]])
}


patient_matrix <- n_genes_p_patients

patient_matrix2 <- patient_matrix # Exclude patients with missing data from clustering

#Get mean gene length
all_lengths <- unlist(lapply(genes_per_patient_names,function(x) length(x)))
media <- mean(all_lengths)

#Obtain optimal clusters
pamk.best <- pamk(patient_matrix2)
kk <- pamk.best$nc ## kk is the optimal number of clusters for the particular optimization.

# Perform hierarchical clustering with the suggested number of clusters
patient_matrix3 <- t(patient_matrix2)
set.seed(2020)
res_hclust <- hclust(jaccard_ind(patient_matrix3),"ward.D2")
res_shc_2 <- shc(t(patient_matrix3), matmet=mfun, linkage="ward.D2", n_sim = 1000,alpha = .08)
res_shc_2$hc_dat$labels <- rownames(t(patient_matrix3))
plot(res_shc_2,alpha=0.8,ci_emp=T,use_labs = TRUE)

# Calculate two 0-1 matrices in order to compare our clustering with the ground truth.
arbol <- cutree(res_hclust,kk)
arbol_splited <- split(names(arbol),arbol)
splited_ground_truth <- split(ground_truth[,1],ground_truth[,2])
arbol_splited_mat <- matrix(0,ncol= nrow(ground_truth),nrow= nrow(ground_truth))
ground_truth_mat <- matrix(0,ncol= nrow(ground_truth),nrow= nrow(ground_truth))
dimnames(arbol_splited_mat) <- list(rownames(ground_truth),rownames(ground_truth))
dimnames(ground_truth_mat) <- list(rownames(ground_truth),rownames(ground_truth))
for(f in 1:nrow(arbol_splited_mat)){
  current_patient_row <- rownames(ground_truth_mat)[f]
  for(g in 1:ncol(arbol_splited_mat)){
    current_patient_col <- colnames(ground_truth_mat)[g]
    cluster_pat_row_ground_truth <- grep(current_patient_row,splited_ground_truth)
    cluster_pat_col_ground_truth <- grep(current_patient_col,splited_ground_truth)
    cluster_pat_row_arbol_splited <- grep(current_patient_row,arbol_splited)
    cluster_pat_col_arbol_splited <- grep(current_patient_col,arbol_splited)
    if(cluster_pat_row_ground_truth==cluster_pat_col_ground_truth){
      ground_truth_mat[f,g] <- 1
      ground_truth_mat[g,f] <- ground_truth_mat[f,g]
    }
    if(cluster_pat_row_arbol_splited==cluster_pat_col_arbol_splited){
      arbol_splited_mat[f,g] <- 1
      arbol_splited_mat[g,f] <- arbol_splited_mat[f,g]
    }
  }
}

sum_matrix <- arbol_splited_mat + ground_truth_mat
tab_sum_matrix <- table(sum_matrix)
zeros <- tab_sum_matrix["0"] # True Negatives
if(is.na(zeros)){
  zeros <- 0
}
twos <- tab_sum_matrix["2"] # True Positives
if(is.na(twos)){
  twos <- 0
}
# Accuracy
accuracy <- (zeros+twos)/sum(tab_sum_matrix,na.rm = T)
final_accuracy_matrix[u,val] <- accuracy
# Matthew's Coefficient
sum_matrix[which(sum_matrix==1,arr.ind=T)] <- 6 #Set all false to value 6.
dif_sum_mat <- sum_matrix - arbol_splited_mat
dif_tab_sum_matrix <- table(dif_sum_mat)
sixes <- dif_tab_sum_matrix["6"] # False Positives
if(is.na(sixes)){
  sixes <- 0
}
fives <- dif_tab_sum_matrix["5"] # False Negatives
if(is.na(fives)){
  fives <- 0
}
tn <- as.double(unname(zeros))
tp <- as.double(unname(twos))
fp <- as.double(unname(fives))
fn <- as.double(unname(sixes))
numerador <- (tp * tn) - (fp * fn)
den_1 <- tp + fp
den_2 <- tp + fn
den_3 <- tn + fp
den_4 <- tn + fn
pro_den <- den_1 * den_2 * den_3 * den_4
denominador <- sqrt(pro_den)
matthews <- numerador/denominador
final_matthews_matrix[u,val] <- matthews
final_kk_used[u,val] <- kk
mean_per_pair[u,val] <- media


genes_interesantes <- lapply(genes_per_patient,names)
genes_interesantes <- lapply(genes_interesantes,function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))

library(gprofiler2)
pathways_res <- lapply(genes_interesantes, function(x) gost(x,'hsapiens'))
pathways_res <- lapply(pathways_res,function(x) x$result)


g1 <- pathways_res[names(which(referencia == 1))]
g1 <- do.call(rbind,g1)
g1 <- names(which(table(g1$term_name)==table(referencia)['1'])) #Ensuring the terms are present in every patient of the group

g1i <- pathways_res[names(which(referencia == 1))]
g1i <- do.call(rbind,g1i)
g1i <- names(which(table(g1i$term_id)==max(table(referencia)['1'])))

genes1 <- genes_interesantes[names(which(referencia == 1))]
genes1 <- do.call(c,genes1)
genes1 <- names(which(table(genes1) == max(table(genes1))))

g2 <- pathways_res[names(which(referencia == 2))]
g2 <- do.call(rbind,g2)
g2 <- names(which(table(g2$term_name)==max(table(referencia)['2'])))

g2i <- pathways_res[names(which(referencia == 2))]
g2i <- do.call(rbind,g2i)
g2i <- names(which(table(g2i$term_id)==max(table(referencia)['2'])))

genes2 <- genes_interesantes[names(which(referencia == 2))]
genes2 <- do.call(c,genes2)
genes2 <- names(which(table(genes2) == max(table(genes2))))

g3 <- pathways_res[names(which(referencia == 3))]
g3 <- do.call(rbind,g3)
g3 <- names(which(table(g3$term_name)==table(referencia)['3']))

g3i <- pathways_res[names(which(referencia == 3))]
g3i <- do.call(rbind,g3i)
g3i <- names(which(table(g3i$term_id)==max(table(referencia)['3'])))

genes3 <- genes_interesantes[names(which(referencia == 3))]
genes3 <- do.call(c,genes3)
genes3 <- names(which(table(genes3) == max(table(genes3))))


lst_all <- list(genes1,genes2,genes3,genes4)
names(lst_all) <- c("Group 1", "Group 2", "Group 3")
saveRDS(lst_all, file="../lst_all")


lst_fun <- list(g1,g2,g3,g4)
names(lst_fun) <- c("Group 1", "Group 2", "Group 3")
lst_id <- list(g1i,g2i,g3i,g4i)
names(lst_id) <- c("Group 1", "Group 2", "Group 3")

saveRDS(lst_fun,"~/Vascular_Disease/genes/function2")
saveRDS(lst_id,"~/Vascular_Disease/genes/id2")
dir.create("~/Vascular_Disease/genes")
capture.output(lst_fun, file="~/Vascular_Disease/genes/function2.txt")
capture.output(lst_id, file="~/Vascular_Disease/genes/id2.txt")


#Venn diagram

ggvenn(
  lst_all, columns=names(lst_all),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

lst_each_2 <- lapply(lst_each,function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'SYMBOL',column = 'ENTREZID')))
ham_f <- lapply(lst_each_2, function(x) hamming_com[x,x])
rownames(ham_f$`Group 1`) <- lapply(rownames(ham_f$`Group 1`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))
rownames(ham_f$`Group 2`) <- lapply(rownames(ham_f$`Group 2`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))
rownames(ham_f$`Group 3`) <- lapply(rownames(ham_f$`Group 3`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))
rownames(ham_f$`Group 4`) <- lapply(rownames(ham_f$`Group 4`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))

colnames(ham_f$`Group 1`) <- lapply(colnames(ham_f$`Group 1`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))
colnames(ham_f$`Group 2`) <- lapply(colnames(ham_f$`Group 2`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))
colnames(ham_f$`Group 3`) <- lapply(colnames(ham_f$`Group 3`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))
colnames(ham_f$`Group 4`) <- lapply(colnames(ham_f$`Group 4`),function(x) unname(mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL')))