library(GSEABase)
library(GSVA)
library(ggplot2)
library(dplyr)
library(viridis)

source("../../dermalNF/bin/dermalNFData.R")
dermals = rna_fpkm_matrix(byIsoform = FALSE)

## step 1 - just combine all
comb = dermals

## step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr <- which(apply(comb,1,function(x) all(x>0.1)))
expr <- setdiff(expr, expr[union(grep("MIR", rownames(comb)[expr]), grep("SNO", rownames(comb)[expr]))])

## step 3, normalize
require(limma)
comb.norm = data.frame(voomWithQualityWeights(comb[expr, ])$E)
comb.norm$Gene = rownames(comb.norm)
comb.norm2 = comb.norm[1:33]
comb.norm2 <- as.matrix(comb.norm2)

colnames(comb.norm2) <- gsub(".", "-", colnames(comb.norm2), fixed = TRUE)
colnames(comb.norm2) <- gsub("X", "", colnames(comb.norm2), fixed = TRUE)

rna.id<-read.table(synGet("syn6156140")@filePath, sep = ",", header = TRUE)
pt.id<-synTableQuery("SELECT 'RNASeq (Cufflinks)', sampleIdentifier FROM syn5556216")@values %>% filter(!is.na(`RNASeq (Cufflinks)`))
colnames(pt.id) <- c("SynapseID", "patientID")
ids <- full_join(rna.id,pt.id)

colnames(comb.norm2)<-ids$patientID[colnames(comb.norm2)==ids$Sample.Id]

## requires gene set collection in wd from MUTSIGDB
## http://software.broadinstitute.org/gsea/downloads.jsp
hallmark.sigs <- getGmt("../../data/h.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"),
                        geneIdType = SymbolIdentifier())

this.file <- "https://github.com/allaway/dNFLandscape/blob/master/analysis/2016-11-22/GSVA_ENCODE_dNF_mini.R"

## hallmarks signature gene set collection GSVA
hallmark.ssGSEA <- gsva(comb.norm2, hallmark.sigs, rnaseq=T, method = "gsva")$es.obs
