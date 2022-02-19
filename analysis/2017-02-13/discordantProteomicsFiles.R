library(synapseClient)
library(dplyr)

new<-read.table(synGet("syn8256971")@filePath, sep = ",", header = TRUE)
old<-read.table(synGet("syn4999547")@filePath, sep = "\t", header = TRUE)

join <- merge(new, old, by.x = "CTF_Sample_Description", by.y = "Patient.ID", all = TRUE)

join <- select(join, CTF_Sample_Description, Sample_No, GSL_Sample_ID, Sample.ID, Sequencing.Library.ID, Proteomics.data, Remark.y)              

names(join) <- c("CTF_Sample_Description", "Sample_No", "New_ID", "Old_ID", "Sequencing_Library_ID", "Proteomics_Data", "Remark")

write.table(join, file = "discordant_proteomics.txt", sep = "\t")

this.file<-"https://raw.githubusercontent.com/allaway/dNFLandscape/200c93be385335f356fb5c9f0a42e882c70bdaf0/analysis/2017-02-13/discordantProteomicsFiles.R"
synStore(File("discordant_proteomics.txt", parentId = "syn7899337"), used = c("syn4999547", "syn8256971"), executed = this.file)
