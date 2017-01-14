library(data.table)
library(synapseClient)
library(dplyr)
library(DEXSeq)

synapseLogin()

TrinityTMM.genes <- synGet('syn7989909')@filePath
Trinity.genes <- fread(file = TrinityTMM.genes, sep = '\t', header = TRUE)

rnaid<-synTableQuery('SELECT RNASeq, sampleIdentifier FROM syn5556216')@values
rnaid<-dplyr::filter(rnaid, !is.na(RNASeq))

rnaid$RNASeq<-lapply(rnaid$RNASeq, function(x) synGet(id=x, downloadFile = FALSE)@fileHandle$fileName)
rnaid$RNASeq<-gsub('accepted_hits_', "", rnaid$RNASeq)
rnaid$RNASeq<-gsub('_featureCounts.txt', "", rnaid$RNASeq)

colnames(Trinity.genes)<-rnaid$sampleIdentifier[match(colnames(Trinity.genes), rnaid$RNASeq)]
rows <- Trinity.genes[,1]
Trinity.genes <- Trinity.genes[,-1]
rownames(Trinity.genes) <- rows$V1
Trinity.genes$means <- rowMeans(Trinity.genes)

