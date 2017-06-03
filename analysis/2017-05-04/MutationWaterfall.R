library(dplyr)
library(tidyr)
library(stringr)
library(synapseClient)
library(GenVisR)

synapseLogin()

mut<-read.table(synGet("syn5839666")@filePath, sep = "\t", header = TRUE, colClasses = "character")
mut$chr <- str_extract(mut$Chromosome, "\\d+")
names(mut)[names(mut)=="Mutation_Type"] <- "Variant_Classification"
names(mut)[names(mut)=="Sample_ID"] <- "Tumor_Sample_Barcode"

soms <- filter(mut, Mutation_Status == "Somatic" & Variant_Classification != "Silent" & Variant_Classification != "Intron")
germs <- filter(mut, Mutation_Status == "Germline" & Variant_Classification != "Silent" & Variant_Classification != "Intron")
set.seed(42)

waterfall(soms, fileType = "MAF", mainDropMut = TRUE, mainRecurCutoff = 0.1, mainPalette = cols)

waterfall(soms %>% filter(Hugo_Symbol %in% cosmic$Gene.Symbol), fileType = "MAF", 
          mainRecurCutoff = 0.0, mainDropMut = TRUE, mainPalette = cols)

waterfall(germs, fileType = "MAF", mainRecurCutoff = 0., mainDropMut = TRUE, mainPalette = cols)

waterfall(germs %>% filter(Hugo_Symbol %in% cosmic$Gene.Symbol), fileType = "MAF", 
          mainRecurCutoff = 0.0, mainDropMut = TRUE, mainPalette = cols)

waterfall(soms, fileType = "MAF", mainDropMut = TRUE, mainRecurCutoff = 0.1, mainPalette = cols, clinData = labs)
