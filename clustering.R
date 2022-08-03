library(sigclust2)

setwd("~/Vascular_Disease")
g.lst <- readRDS("~/Vascular_Disease/sel_genes_per_patient")

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


datalist <- list()

for (i in 1:length(g.lst)) {
  v <- c()
  for (j in 1:length(g.lst)) {
    v <- c(v, jaccard(g.lst[[i]],g.lst[[j]]))
  }
  datalist[[i]] <- v
}  

datadf <- do.call(rbind,datalist)
row.names(datadf) <- names(g.lst)
colnames(datadf) <- names(g.lst)
diag(datadf) <- 0

shc_result <- shc(datadf, metric="euclidian", linkage="ward.D2")

png(file="plots/hc/hc_pre_multi.png", width =465, height = 225, units = "mm", res=300)

plot(shc_result, hang=.1)

dev.off()
