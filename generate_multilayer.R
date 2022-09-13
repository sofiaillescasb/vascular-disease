library(CmmD)

# Community trajectory matrix and hamming distance matrix were originally generated with these parameters
# Create a vector with the paths where Molti's Output files are saved.
structures_12 <- paste0("data/Molti_Output/",seq(0.5,12,0.5),".csv")
# Detect community trajectories and tree distances between each gene.
matrices <- CmmD_from_community_structures(nodelist = NULL, community_structures = structures_12, resolution_start = 0.5,resolution_end = 12,interval = 0.5,distmethod = "hamming",threads = 7)
matrices$hamming_distance_matrix = as.matrix(matrices$distance_matrix) * 24 # This transformation is needed because parallel dist is weighted.
# 24 = length(seq(0.5,12,0.5)) -> number of resolution values analyzed