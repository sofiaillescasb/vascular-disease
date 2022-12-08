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

colnames(all_samples)[-1] <- c("HUVECS", "VM015",  "VM024",  "VM038",  "VM040",  "VM042",  "VM043",  "VM048",  "VM053",  "VM054",  "VM055", 
                          "VM056",  "VM060",  "VM064",  "VM066",  "VM068",  "VM071",  "VM072",  "VM073",  "VM081",  "VM082",  "VM083", 
                           "VM085",  "VM089",  "VM090",  "VM092",  "VM093",  "VM099",  "VM103",  "VM108",  "VM110",  "VM111",  "VM113", 
                           "VM119",  "VM124", "VM125",  "VM127")

write.csv(all_samples, "../counts/assay.csv", row.names = FALSE, quote = FALSE)

#Importing counts data
counts <- read.csv("../counts/assay.csv", header=TRUE, row.names=1)
list(counts)
head(rownames(counts))
#Importing sample metadata
meta <- read.csv("../se/MatrzMeta.csv", header=TRUE, row.names = 2)
head(meta)

#Creating SummarizedExperiment object
se <- SummarizedExperiment(assays=list(counts=counts), colData=meta)

# Data exploration
metadata <- as.data.frame(colData(se))
head(metadata)

saveRDS(se, file="../se/se_org")
