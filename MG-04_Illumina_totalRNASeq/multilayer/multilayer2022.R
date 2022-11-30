library(CmmD)

setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq")
#done with the marenostrum
# Create a vector with the paths where Molti's Output files are saved. 
structures_12 <- paste0("out_molti/",seq(0,50,1),".csv")
# Detect community trajectories and tree distances between each gene. 
curie_to_12_full <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0,resolution_end = 50,interval = 1,distmethod = "hamming",threads = 7)
curie_to_12_full$hamming_distance_matrix = as.matrix(curie_to_12_full$distance_matrix) * 50 # This transformation is needed because parallel dist is weighted.
write.csv(curie_to_12_full$hamming_distance_matrix, file="multilayer/hamming_2022.csv")
