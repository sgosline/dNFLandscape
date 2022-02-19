library(dplyr)
library(tidyr)
library(stringr)
library(synapseClient)
library(GenVisR)

this.file = "https://github.com/allaway/dNFLandscape/blob/master/analysis/2017-05-04/MutationWaterfall.R"

synapseLogin()
cosmic <- read.table(synGet("syn9754614")@filePath, sep = "\t", header = TRUE, colClasses = "character")
mut<-read.table(synGet("syn5839666")@filePath, sep = "\t", header = TRUE, colClasses = "character", stringsAsFactors = TRUE)

mut$chr <- str_extract(mut$Chromosome, "\\d")
names(mut)[names(mut)=="Mutation_Type"] <- "Variant_Classification"
names(mut)[names(mut)=="Sample_ID"] <- "Tumor_Sample_Barcode"

mut$cosmic[mut$Hugo_Symbol %in% cosmic$Gene.Symbol] <- "T"
mut$cosmic[!mut$Hugo_Symbol %in% cosmic$Gene.Symbol] <- "F"

soms <- filter(mut, Mutation_Status == "Somatic" & Variant_Classification != "Silent" & Variant_Classification != "Intron" & Variant_Classification != "In_Frame_Ins" & Variant_Classification != "In_Frame_Del")

germs <- filter(mut, Mutation_Status == "Germline" & Variant_Classification != "Silent" & Variant_Classification != "Intron" & Variant_Classification != "In_Frame_Ins" & Variant_Classification != "In_Frame_Del")

set.seed(42)

cols <- c("#17B890", "#FF6F59", "#D7263D", "black", "black", "black", "black", "#A7BBEC", "#FFF07C", "black", "black", "black", "black", "black", "black", "black", "black", "black")


###All somatic mutations with over 10% occurrence rate
svg(file = "Soms10perc.svg", height = 11, width = 10)
waterfall(soms, fileType = "MAF", mainDropMut = TRUE, mainRecurCutoff = 0.1, mainPalette = cols)
dev.off()
synStore(File("Soms10perc.svg", parentId = "syn10968238"), used = c("syn9754614", "syn5839666"), executed = this.file)

###COSMIC listed genes with somatic mutations, no cutoff
svg(file = "SomsCOSMICFilt.svg", height = 11, width = 10)
waterfall(soms[soms$Hugo_Symbol %in% cosmic$Gene.Symbol,], fileType = "MAF", 
          mainRecurCutoff = 0.0, mainDropMut = TRUE, mainPalette = cols)
dev.off()
synStore(File("SomsCOSMICFilt.svg", parentId = "syn10968238"), used = c("syn9754614", "syn5839666"), executed = this.file)


###All germline mutations, no cutoff
svg(file = "Germs80perc.svg", height = 11, width = 10)
waterfall(germs, fileType = "MAF", mainRecurCutoff = 0.7, mainDropMut = TRUE, mainPalette = cols)
dev.off()
synStore(File("Germs80perc.svg", parentId = "syn10968238"), used = c("syn9754614", "syn5839666"), executed = this.file)

###COSMIC listed genes with germline mutations, no cutoff
svg(file = "GermsCOSMICFilt.svg", height = 11, width = 10)
waterfall(germs[germs$Hugo_Symbol %in% cosmic$Gene.Symbol,], fileType = "MAF", 
          mainRecurCutoff = 0.0, mainDropMut = TRUE, mainPalette = cols)
dev.off()
synStore(File("GermsCOSMICFilt.svg", parentId = "syn10968238"), used = c("syn9754614", "syn5839666"), executed = this.file)


soms <- filter(mut, Mutation_Status == "Germline") %>% filter(Hugo_Symbol %in% c("CDC27", "CREBBP"))
svg(file = "CREBBPCDC27perc.svg", height = 11, width = 10)
waterfall(soms, fileType = "MAF", mainDropMut = TRUE, mainPalette = cols)
dev.off()

synStore(File("CREBBPCDC27perc.svg", parentId = "syn10968238"), used = c("syn9754614", "syn5839666"), executed = this.file)
