library(data.table)
library(SummarizedExperiment)

setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/counts")
tmp = list.files(pattern = "out.tab")
tmp
lst <- lapply(tmp, read.table)
lst <- lapply(lst, "[", -c(1:4),)
files_matrix <- rbind(lst)
files_matrix <- lapply(files_matrix, setNames, nm = c("SAMPLE", "UNSTRANDED", "FORWARD", "REVERSE"))
files_matrix[[1]][1:10,]


#Adding all gene counts in the reverse column 
reverse <- sum(files_matrix[[length(files_matrix)]]$REVERSE[-(1:5)])

#doing the same for forward column to see which has a higher quantity

forward <- sum(files_matrix[[length(files_matrix)]]$FORWARD[-(1:5)])

#Use column with higher counts

if (reverse > forward) {
  higher_counts <- lapply(files_matrix, "[[", "REVERSE")
} else if (forward > reverse) {
  higher_counts <- lapply(files_matrix, "[[", "FORWARD")
} else {
  higher_counts <- lapply(files_matrix, "[[", "UNSTRANDED")
}

gene_names <- files_matrix[[1]]$SAMPLE
all_samples <- data.frame(GENE = gene_names)
samples <- c("HUVECs_2", "HUVECs_3", "HUVECs_4", "HDLECs_1", "HDLECs_2", "HDLECs_3")

for (i in 2:(length(files_matrix)+1)) {
  all_samples[i] = higher_counts[i-1] 
}


colnames(all_samples)[-1] <- samples
head(all_samples)
write.csv(all_samples, "../preprocess/assay.csv", row.names = FALSE, quote = FALSE)

patients <- read.csv("../../counts/assay.csv", header=TRUE)
setdiff(patients[1],all_samples[1]) 

counts <- read.csv("../preprocess/assay.csv", header=TRUE, row.names = 1)
#Importing sample metadata
meta <- read.csv("../preprocess/MatrzMeta.csv", header=TRUE, row.names = 2)
head(meta)

#Creating SummarizedExperiment object
se <- SummarizedExperiment(assays=list(counts=counts), colData=meta)
saveRDS(se, file="../preprocess/se_contrl")


