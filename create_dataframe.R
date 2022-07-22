library(data.table)
library(SummarizedExperiment)

setwd("~/Vascular_Disease/Reads_Per_Gene")
tmp = list.files(pattern = "out.tab")
tmp
lst <- lapply(tmp, read.table)
lst <- lapply(lst, "[", -c(1:4),)
files_matrix <- rbind(lst)
files_matrix <- lapply(files_matrix, setNames, nm = c("SAMPLE", "UNSTRANDED", "FORWARD", "REVERSE"))
files_matrix[[1]][1:10,]


#Adding all gene counts in the reverse column 
reverse <- sum(files_matrix[[37]]$REVERSE[-(1:5)])

#doing the same for forward column to see which has a higher quantity

forward <- sum(files_matrix[[37]]$FORWARD[-(1:5)])

#Reverse has a higher count, use 4th column

reverse_counts <- lapply(files_matrix, "[[", "REVERSE")
gene_names <- files_matrix[[1]]$SAMPLE
all_samples <- data.frame(GENE = gene_names)
samples <- read.table("/home/bscuser/Vascular_Disease/list.txt")
samples

for (i in 2:38) {
  all_samples[i] = reverse_counts[i-1] 
}

colnames(all_samples) <- unlist(samples)
all_samples
write.csv(all_samples, "assay.csv", row.names = FALSE, quote = FALSE)
