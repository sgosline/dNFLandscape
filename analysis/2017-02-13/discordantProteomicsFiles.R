library(synapseClient)
library(dplyr)

new<-read.table(synGet("syn8256971")@filePath, sep = ",", header = TRUE)
old<-read.table(synGet("syn4999547")@filePath, sep = "\t", header = TRUE)

join <- merge(new, old, by.x = "CTF_Sample_Description", by.y = "Patient.ID", all = TRUE)

join <- select(join, CTF_Sample_Description, Sample_No, GSL_Sample_ID, Sample.ID, Sequencing.Library.ID, Proteomics.data, Remark.y)              

names(join) <- c("CTF_Sample_Description", "Sample_No", "New_ID", "Old_ID", "Sequencing_Library_ID", "Proteomics_Data", "Remark")

synStore()