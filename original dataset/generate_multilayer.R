setwd('/home/bscuser/Vascular_Disease/Medulloblastoma-main')

library(CmmD)

# The following script is intended to be run at clusters running SLURM and Greasy. It is programed to be run from bash in the following way:
# RScript Randomized_curie.R 1 50 
# With such command, the script will perform randomizations 1 to 50 and save 2 files, the accuracy matrix and the number of clusters suggested for each theta-lambda pair.

set.seed(NULL) # Restart all seeds (Pure randomizing)

# Load Dependencies
library(pvclust)
library(fpc)
library(readr)
library(knitr)
library("AnnotationDbi")
library("igraph")
library("stringr")
library("e1071")

# Define Jaccard index function to use with hclust
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
  return(as.dist(1-res))
}

# Other useful function
unsplit_to_data_frame <- function(x){
  unlisted <- unlist(x,use.names=T)
  col_1 <- unlist(x,use.names=F)
  col_2 <- substr(names(unlisted),start=1,stop=4)
  tabla <- as.data.frame(matrix(nrow=length(unlisted),ncol=2))
  tabla[,1] <- col_1
  tabla[,2] <- col_2
  return(tabla)
}

# Define CmmD Package functions (Not installed at our supercomputer)
CmmD <- function (nodelist = NULL, input_layers, resolution_start, resolution_end,
                  interval, destfile_community_analysis)
{
  require("AnnotationDbi")
  require("igraph")
  require("stringr")
  require("e1071")
  if (length(input_layers) < 2) {
    stop("ERROR: Input_layers argument must be a list of at least 2 network files")
  }
  if (class(resolution_end) != "numeric") {
    stop("ERROR: Resolution parameter must be a number")
  }
  if (class(resolution_start) != "numeric") {
    stop("ERROR: Resolution parameter must be a number")
  }
  if (class(interval) != "numeric") {
    stop("ERROR: Interval value must be a number")
  }
  if (class(destfile_community_analysis) != "character") {
    stop("ERROR: destfile_community_analysis expects a character string")
  }
  layers <- paste0(input_layers, collapse = " ")
  message(paste0("Resolution parameter starts at: ", resolution_start))
  message(paste0("Resolution parameter ends at: ", resolution_end))
  resolution_interval <- seq(from = resolution_start, to = resolution_end,
                             by = interval)
  desfile_vector <- paste0(destfile_community_analysis, resolution_interval,
                           ".csv")
  message(paste0("Starting community analysis."))
  start_time <- Sys.time()
  for (i in 1:length(resolution_interval)) {
    current_resolution <- resolution_interval[i]
    current_destfile <- desfile_vector[i]
    current_layers <- layers
    message(paste0("Resolution parameter: ", current_resolution))
    message(Sys.time())
    system_order <- paste("molti-console", "-o", current_destfile,
                          "-p", current_resolution, layers)
    system(system_order)
  }
  message(paste0("Reading MolTi output files. Calculating Gene/Community matrix"))
  output_files <- list.files(destfile_community_analysis)
  to_be_forgotten <- grep("_", output_files)
  output_files <- output_files[-to_be_forgotten]
  alllists <- list()
  for (i in 1:length(output_files)) {
    red <- readLines(paste0(destfile_community_analysis,
                            output_files[i]))
    cluster_ids <- grep("Cluster", red)
    lista <- list()
    for (j in 1:length(cluster_ids)) {
      st <- cluster_ids[j]
      if (j == length(cluster_ids)) {
        en <- length(red)
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-length(current_cluster)]
      }
      else {
        en <- cluster_ids[j + 1]
        current_cluster <- red[st:en]
        current_cluster2 <- current_cluster[-c(length(current_cluster),
                                               length(current_cluster) - 1)]
      }
      lista[[j]] <- current_cluster2[2:length(current_cluster2)]
      names(lista)[j] <- paste0("Cluster_", j)
    }
    kaz <- output_files[i]
    assign(paste0("com_", kaz), value = lista)
    alllists[[i]] <- lista
  }
  names(alllists) <- output_files
  tamano_alllists <- length(alllists)
  allgenes <- unique(unlist(alllists))
  if (length(nodelist) > 0) {
    inter_nodes <- intersect(allgenes, nodelist)
    allgenes <- inter_nodes
  }
  print(paste0("Files red. Calculating Gene/Community matrix"))
  res_matrix <- matrix(ncol = tamano_alllists + 1, nrow = length(allgenes))
  rownames(res_matrix) <- allgenes
  colnames(res_matrix) <- c(output_files, "Pattern")
  for (i in 1:length(allgenes)) {
    gen <- rownames(res_matrix)[i]
    for (j in 1:tamano_alllists) {
      searched <- unlist(lapply(alllists[[j]], function(x) gen %in%
                                  x))
      comunidad <- unname(which(searched == TRUE))
      res_matrix[i, j] <- comunidad
    }
    res_matrix[i, "Pattern"] <- paste0(res_matrix[i, 1:(ncol(res_matrix) -
                                                          1)], collapse = "_")
    percentage <- round((i/length(allgenes)), digits = 4) *
      100
    porcentajes <- seq(from = 0, to = 100, by = 5)
    progresos <- paste0("Progress: ", porcentajes, "%")
    porcentajes[1] <- 1
    if ((percentage %in% porcentajes) == TRUE) {
      cual_percentage <- which(porcentajes == percentage)
      to_post <- paste0("Progress: ", percentage, "%")
      message(to_post)
    }
  }
  message(paste0("Gene/Community matrix calculated, calculating Hamming distances for all gene pairs. This process may take a while: It takes about 14 min with an Intel Xeon E-2124 processor"))
  genes_same_communities <- split(rownames(res_matrix), res_matrix[,
                                                                   "Pattern"])
  final_res_matrix_length <- ncol(res_matrix) - 1
  distance_matrix <- hamming.distance(res_matrix[, 1:final_res_matrix_length])
  final_output <- list(res_matrix[, 1:final_res_matrix_length],
                       genes_same_communities, distance_matrix)
  names(final_output) <- c("gene_community_matrix", "l_constant",
                           "hamming_distance_matrix")
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  message(paste0("Run Time: ", diff_time))
  return(final_output)
}
nets <- list.files("/home/bscuser/Vascular_Disease/gene_multilayer_network-master/networks/")
m.n.2022 <- CmmD(nodelist = NULL,input_layers=nets,resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = "hamming",threads = 7,destfile_community_analysis = "/home/bscuser/Vascular_Disease/Medulloblastoma-main/data/Molti_Output")
