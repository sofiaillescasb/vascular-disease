library(knitr)
library(RColorBrewer)
library(SummarizedExperiment)
library(geneplotter)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(edgeR)
library(corpcor)  
library(M3C)
library(calibrate)
library(EnhancedVolcano) 
library(sva) 
library(pvclust) # Hierarchical clustering with p-values
library(dendextend) # fancy plotting dendograms
library(caret)
library(psych)
library(sigclust2)

setwd("~/Vascular_Disease")

#Importing counts data
counts <- read.csv("counts/assay.csv", header=TRUE, row.names=1)
list(counts)
head(rownames(counts))
#Importing sample metadata
meta <- read.csv("se/MatrzMeta.csv", header=TRUE, row.names = 2)
head(meta)

#Creating SummarizedExperiment object
se <- SummarizedExperiment(assays=list(counts=counts), colData=meta)

# Data exploration
metadata <- as.data.frame(colData(se))
fplot <- ggplot(as.data.frame(metadata))
head(metadata)


#Exploring metadata
# 
fplot + geom_bar(aes(Mutant.gene, fill = Variant)) + xlab("Mutation") + ylab("Number of Samples") +
  ggtitle("Sample mutation distribution") + scale_fill_brewer(palette = "Set3")

fplot + geom_bar(aes(Summary.clinic, fill = Mutant.gene)) + xlab("Summary clinic") + ylab("Number of Samples") +
  ggtitle("Sample summary clinic distribution") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_fill_brewer(palette = "Set3")

fplot + geom_bar(aes(Summary.pathological.marker, fill = Mutant.gene)) + xlab("Pathological marker") +
  ylab("Number of Samples") + ggtitle("Sample pathological marker distribution") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_brewer(palette = "Set3")

fplot + geom_bar(aes(Summary.pathological.marker, fill = Variant)) + xlab("Pathological marker") +
  ylab("Number of Samples") + ggtitle("Sample pathological marker distribution") +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + scale_fill_brewer(palette = "Set3")

# Quality assessment and normalization

## Preprocessing data
dge <- DGEList(counts = assays(se)$counts, genes = as.data.frame(mcols(se)), group = se$Case.Control)

# Computing counts per million
assays(se)$logCPM <- cpm(dge, log=TRUE, prior.count=0.5)
assays(se)$logCPM[5:10, 5:10]

dge$counts[5:10, 5:10]
dge$samples$lib.size

## Examining the sequencing depth
ord <- order(dge$samples$lib.size)
dge$samples$lib.size[ord]


# Barplot patients/controls
barplot(dge$samples$lib.size[ord]/1e+06, las = 1, main="Barplot - library size", ylab = "Million of counts", xlab = "Samples", col = c("blue","red")[(se$Case.Control[ord] == "case") +1], border = NA)
legend("topleft", c("cases","control"), fill = c("red","blue"), inset = 0.01)

# QQplot
sampledepth <- round(dge$sample$lib.size / 1e6, digits=1)
qqnorm(sampledepth)
qqline(sampledepth)

# Distribution of expression levels among genes
avgexp <- rowMeans(assays(se)$logCPM)
hist(avgexp, xlab="log2 CPM", main="Distribution of the genes average expression", las=1)
abline(v=1, col="red", lwd=2)


# Filtering lowly expressed genes
nsamples <- length(se$Sample.ID)
sample_cutoff <- 0.2
nsamples_cutoff <- sample_cutoff*nsamples
nsamples_cutoff

logcpm_cutoff <- 1
mask <- rowSums(assays(se)$logCPM <= logcpm_cutoff) >= nsamples_cutoff

dim(se)
se.filt <- se[!mask, ]
dim(se.filt)
dge.filt <- dge[!mask, ]
dim(dge.filt)
kept_genes2b <- dim(se.filt)[1]

# Between sample normalization
dge.filt.norm <- calcNormFactors(dge.filt)
assays(se.filt)$logCPM.norm <- cpm(dge.filt.norm, log = TRUE, prior.count = 3, normalized.lib.sizes = TRUE)

# Plots: in-sample normalized only, filtered, filtered and normalized between samples

multidensity(as.list(as.data.frame(assays(se[, se$Case.Control == "case"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)

multidensity(as.list(as.data.frame(assays(se.filt[, se.filt$Case.Control == "case"])$logCPM.norm)),
             xlab="log 2 CPM", legend=NULL, main="Patient samples", las=1)


# PCA based on normalized counts(normalized between samples)
pca_norm <- prcomp(as.data.frame(t(assays(se.filt)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm, var.axes=FALSE, groups=se.filt$Case.Control, labels = se.filt$Sample.ID) + ggtitle("PCA based on normalized data and colored by type")

# Filtering the three outliers

dge.filt.norm_2 <- dge.filt.norm[,-c(3,5,28)]
se.filt_2 <- se.filt[,-c(3,5,28)]
metadata_2 <- metadata[-c(3,5,28),]

# PCA based on normalized counts(normalized between samples) excluding outlying samples

pca_norm_2 <- prcomp(as.data.frame(t(assays(se.filt_2)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_2, var.axes=FALSE, groups=se.filt_2$Case.Control, labels = se.filt_2$Sample.ID) + ggtitle("PCA based on normalized data and colored by type")

#t-SNE
##t-SNE based on normalized counts
tsne_by_perplex <- function(s=se.filt, metadata=metadata) {
  
  nsamples <- length(s$Sample.ID)
  max_perplex <- trunc(((nsamples - 1)/3)-1)
  optimal_perplex = round(nsamples^(1/2))
  perplex1 = round(optimal_perplex/2)
  perplex3 = mean(c(optimal_perplex,max_perplex))
  
  s$simplified_type <- s$Case.Control
  s$simplified_type <- gsub("control","C",s$simplified_type)
  s$simplified_type <- gsub("case","",s$simplified_type)
  
  if(perplex1 <= max_perplex){
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Case.Control, perplex=perplex1, seed=TRUE) + ggtitle("tsne based on normalized counts perplex1"))
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Mutant.gene, perplex=perplex1, seed=TRUE, text=as.character(s$simplified_type)) + ggtitle("tsne based on normalized counts perplex1"))
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Variant, perplex=perplex1, seed=TRUE) + ggtitle("tsne based on normalized counts perplex1"))
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Summary.clinic, perplex=perplex1, seed=TRUE) + ggtitle("tsne based on normalized counts perplex1"))
  
  
  if(optimal_perplex <= max_perplex){
    print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Case.Control, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne based on normalized counts perplex2"))
    print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Mutant.gene, perplex=optimal_perplex, seed=TRUE, text=as.character(s$simplified_type)) + ggtitle("tsne based on normalized counts perplex2"))
    print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Variant, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne based on normalized counts perplex2"))
    print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Summary.clinic, perplex=optimal_perplex, seed=TRUE) + ggtitle("tsne based on normalized counts perplex2"))
    
    
    if(perplex3 <= max_perplex){
      print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Case.Control, perplex=perplex3, seed=TRUE) + ggtitle("tsne based on normalized counts perplex3"))
      print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Mutant.gene, perplex=perplex3, seed=TRUE, text=as.character(s$simplified_type)) + ggtitle("tsne based on normalized counts perplex3"))
      print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Summary.clinic, perplex=perplex3, seed=TRUE) + ggtitle("tsne based on normalized counts perplex3"))
      print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Variant, perplex=perplex3, seed=TRUE) + ggtitle("tsne based on normalized counts perplex3"))
    }
  }
}

if(max_perplex %in% c(perplex1,optimal_perplex,perplex3)){
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Case.Control, perplex=max_perplex, seed=TRUE)  + ggtitle("tsne based on normalized counts max_perplex"))
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Mutant.gene, perplex=max_perplex, seed=TRUE, text=as.character(s$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Summary.clinic, perplex=max_perplex, seed=TRUE) + ggtitle("tsne based on normalized counts perplex1"))
  print(tsne(as.data.frame(assays(s)$logCPM.norm), labels=s$Variant, perplex=max_perplex, seed=TRUE, text=as.character(s$simplified_type)) + ggtitle("tsne based on logCPM perplex3"))
}
}
tsne_by_perplex()

## t-SNE based on normalized counts excluding outlying samples

tsne_by_perplex(se.filt_2, metadata_2)


# spectrum clustering based on normalized counts

# Hierarchical clustering based on filtered and normalized counts using spearman correlation

## Hierarchical clustering

plot_hierarchichal_clustering <- function(data,info, save_clusters=FALSE, key_word="unkown", n_boot=100){
  
  # pvclust
  pvclust_result <- pvclust(as.matrix(data), method.hclust = "average", method.dist = "correlation", nboot=n_boot, parallel=TRUE, iseed = 1)
  
  # Basic
  par(mfrow=c(1, 1))
  plot(pvclust_result)
  clusters <- pvpick(pvclust_result) # Returns the clusters and GSMs in it.
  pvrect(pvclust_result)

  # # Coloring by mutation
  colors = brewer.pal(length(unique(metadata$Mutant.gene)),"Set3")
  colors = colors[as.factor(se.filt$Mutant.gene)]
  dend <- as.dendrogram(pvclust_result)
  colors = colors[order.dendrogram(dend)]
  labels_colors(dend)<-rep(0,nsamples) # Changing the color of the text labels to white
  dend %>% set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", 0.6) %>%  # node point size
    # set("leaves_col", colors[as.factor(se.filt$series)]) %>%
    set("leaves_col", colors) %>%
    # set("leaves_col", colors[as.integer(factor(se.filt$type))]) %>% # node point color
    plot(main = paste("Hierarchichal clustering",info,sep=" "))
  legend("topright",legend = unique(se.filt$Mutant.gene[order.dendrogram(dend)]), fill = unique(colors), cex=0.6, box.lty=0)
  pvrect(pvclust_result)
  
  # # Coloring by variant
  colors = brewer.pal(length(unique(metadata$Variant)),"Set3")
  colors = colors[as.factor(se.filt$Variant)]
  colors = colors[order.dendrogram(dend)]
  labels_colors(dend)<-rep(0,nsamples) # Changing the color of the text labels to white
  dend %>% set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", 0.6) %>%  # node point size
    # set("leaves_col", colors[as.factor(se.filt$series)]) %>%
    set("leaves_col", colors) %>%
    # set("leaves_col", colors[as.integer(factor(se.filt$type))]) %>% # node point color
    plot(main = paste("Hierarchichal clustering",info,sep=" "))
  legend("topright",legend = unique(se.filt$Variant[order.dendrogram(dend)]), fill = unique(colors), cex=0.6, box.lty=0)
  pvrect(pvclust_result)
  
  # # Coloring by summary clinic
  colors = c(brewer.pal(length(unique(metadata$Summary.clinic)[1:8]),"Set3"), brewer.pal(length(unique(metadata$Summary.clinic)[9:16]),"Set2"))
  colors = colors[as.factor(se.filt$Summary.clinic)]
  colors = colors[order.dendrogram(dend)]
  labels_colors(dend)<-rep(0,nsamples) # Changing the color of the text labels to white
  dend %>% set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", 0.6) %>%  # node point size
    # set("leaves_col", colors[as.factor(se.filt$series)]) %>%
    set("leaves_col", colors) %>%
    # set("leaves_col", colors[as.integer(factor(se.filt$type))]) %>% # node point color
    plot(main = paste("Hierarchichal clustering",info,sep=" "))
  legend("topright",legend = unique(se.filt$Summary.clinic[order.dendrogram(dend)]), fill = unique(colors), cex=0.6, box.lty=0)
  pvrect(pvclust_result)
  
  # # Coloring by pathological marker
  colors = c(brewer.pal(length(unique(metadata$Summary.pathological.marker)[1:7]),"Set3"), brewer.pal(length(unique(metadata$Summary.clinic)[8:14]),"Set2"))
  colors = colors[as.factor(se.filt$Summary.pathological.marker)]
  colors = colors[order.dendrogram(dend)]
  labels_colors(dend)<-rep(0,nsamples) # Changing the color of the text labels to white
  dend %>% set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", 0.6) %>%  # node point size
    # set("leaves_col", colors[as.factor(se.filt$series)]) %>%
    set("leaves_col", colors) %>%
    # set("leaves_col", colors[as.integer(factor(se.filt$type))]) %>% # node point color
    plot(main = paste("Hierarchichal clustering",info,sep=" "))
  legend("topright",legend = unique(se.filt$Summary.pathological.marker[order.dendrogram(dend)]), fill = unique(colors), cex=0.6, box.lty=0)
  pvrect(pvclust_result)
  
  return(clusters)
}

shc1 <- shc(as.matrix(as.data.frame(t(assays(se.filt)$logCPM.norm))), linkage="ward.D2", n_sim = 100)
plot(shc1,alpha=0.5,ci_emp=T,use_labs = TRUE)

logCPM <- cpm(dge.filt.norm, log=TRUE, prior.count=3)
#normalized_clusters <- plot_hierarchical_clustering(logCPM,'based on normalized counts', n_boot=1000)

## Hierarchical clustering no outliers

logCPM_2 <- cpm(dge.filt.norm_2, log=TRUE, prior.count=3)
normalized_clusters <- plot_hierarchical_clustering(logCPM_2,'based on normalized counts', n_boot=1000)


## Finding highly variable genes

v <- apply(assay(se.filt), 1, var)
df <- as.data.frame(v)
colnames(df) <- 'Variance'


sel <- order(v, decreasing=TRUE)[1:1100]
se.sel <- se.filt[sel]
dge.sel <- DGEList(counts = assays(se.sel)$counts, genes = as.data.frame(mcols(se.sel)), group = se.sel$Case.Control)


#PCA based only on genes with highest variance

pca_norm_sel <- prcomp(as.data.frame(t(assays(se.sel)$logCPM.norm)) , center=TRUE, scale=TRUE)
ggbiplot(pca_norm_sel, var.axes=FALSE, groups=se.sel$Case.Control, labels = se.sel$Sample.ID) + ggtitle("PCA genes with higher variability")

#t-SNE based on genes with highest variance
tsne_by_perplex(se.sel)

#Hierarchical clustering based on genes with highest variance

logCPM_sel <- cpm(dge.sel, log=TRUE, prior.count=3)
normalized_clusters <- plot_hierarchichal_clustering(logCPM_sel,'based on normalized counts for genes with highest variance', n_boot=1000)

counts.sel <- as.data.frame(dge.sel$counts)
row.names(counts.sel) <- rownames(dge.sel$
                                    counts)
write.csv(counts.sel, "~/Vascular_Disease/counts_sel_overmean.csv")

##Finding correlated genes, creating correlation matrix

c <- corr.test(t(assay(se.sel)), method="spearman")
corr.df <- as.data.frame(c$r)
corr.df.p <- as.data.frame(c$p)
corr.df.t <- as.data.frame(c$t)
write.csv(corr.df, "~/Vascular_Disease/correlation_matrix.csv")
write.csv(corr.df.p, "~/Vascular_Disease/correlation_p.csv")
write.csv(corr.df.t, "~/Vascular_Disease/correlation_t.csv")

