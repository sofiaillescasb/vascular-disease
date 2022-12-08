library(AnnotationDbi)
library(org.Hs.eg.db)

same_traj_2 <- lapply(same_traj,function(x) mapIds(x = org.Hs.eg.db,keys =x,keytype = 'ENTREZID',column = 'SYMBOL'))

selg <- readRDS("~/Vascular_Disease/sel_genes_per_patient")
lst <- lapply(selg, function(x) intersect(x,same_traj_2))
lst <- unlist(lst)

corr.comp$membership[n.sel]