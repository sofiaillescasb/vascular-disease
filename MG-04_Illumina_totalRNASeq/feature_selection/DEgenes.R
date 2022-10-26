setwd("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/feature_selection")

library(ggvenn)

se.filt <- readRDS("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/se_filt")
genes_lst <- readRDS("~/Vascular_Disease/MG-04_Illumina_totalRNASeq/preprocess/DEgenes")

#Comparing all controls to all patients
patients <- assays(se.filt)$logCPM.norm[genes_lst$`all patients vs all controls`,c(2:37)]
controls <- assays(se.filt)$logCPM.norm[genes_lst$`all patients vs all controls`,-c(2:37)]

c_mean <- as.matrix(apply(controls, 1, mean))
colnames(c_mean) <- "Mean Expression"

fc <- apply(patients, 2, function(x) abs(x-c_mean))
rownames(fc) <- rownames(patients)

all_per_patient <- apply(fc, 2, function(x) x[x>1])
genes_per_patient <- lapply(all_per_patient, function(x) names(x))


#Comparing venous malformation patients to venous controls and lymphatic patients to lymphatic controls
patients_v <- assays(se.filt[genes_lst$venous,se.filt$Summary.clinic=="venous malformation"])$logCPM.norm
patients_l <- assays(se.filt[genes_lst$lymphatic,se.filt$Summary.clinic=="lymphatic malformation"])$logCPM.norm

controls_v <- assays(se.filt)$logCPM.norm[genes_lst$venous,c("HUVECS", "HUVECs_2", "HUVECs_3", "HUVECs_4")]
controls_l <- assays(se.filt)$logCPM.norm[genes_lst$lymphatic,c("HDLECs_1", "HDLECs_2", "HDLECs_3")]

mean_v <- as.matrix(apply(controls_v, 1, mean))
mean_l <- as.matrix(apply(controls_l, 1, mean))

fc_v <- apply(patients_v, 2, function(x) abs(x-mean_v))
rownames(fc_v) <- rownames(patients_v)
fc_l <- apply(patients_l, 2, function(x) abs(x-mean_l))
rownames(fc_l) <- rownames(patients_l)

v_per_patient <- apply(fc_v, 2, function(x) x[x>1])
genes_per_patient_v <- lapply(v_per_patient, function(x) names(x))
l_per_patient <- apply(fc_l, 2, function(x) x[x>1])
genes_per_patient_l <- rownames(l_per_patient)

#Trying to compare the values of DEgenes in venous control with the DE values of all the patients to try to find differences with non venous patients 
not_in_common <- setdiff(genes_lst$venous, genes_lst$`all patients vs all controls`)
comp_patients <- assays(se.filt)$logCPM.norm[genes_lst$venous[genes_lst$venous%in%not_in_common],c(2:37)]
fc_comp_v <- apply(comp_patients, 2, function(x) abs(x-mean_v[rownames(mean_v)[rownames(mean_v)%in%not_in_common],]))
sig_comp_v <- apply(fc_comp_v, 2, function(x) x[x>1])
