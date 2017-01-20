library(data.table)
library(synapseClient)
library(dplyr)
library(dexus)

synapseLogin()

TrinityTMM.trans <- synGet('syn8034127')@filePath
Trinity.trans <- fread(file = TrinityTMM.trans, sep = '\t', header = TRUE)

rnaid<-synTableQuery('SELECT RNASeq, sampleIdentifier FROM syn5556216')@values
rnaid<-dplyr::filter(rnaid, !is.na(RNASeq))

rnaid$RNASeq<-lapply(rnaid$RNASeq, function(x) synGet(id=x, downloadFile = FALSE)@fileHandle$fileName)
rnaid$RNASeq<-gsub('accepted_hits_', "", rnaid$RNASeq)
rnaid$RNASeq<-gsub('_featureCounts.txt', "", rnaid$RNASeq)

colnames(Trinity.trans)<-rnaid$sampleIdentifier[match(colnames(Trinity.trans), rnaid$RNASeq)]
rows <- Trinity.trans[,1]
Trinity.trans <- Trinity.trans[,-1]
Trinity.trans2 <- as.matrix(Trinity.trans)
rownames(Trinity.trans2) <- rows$V1

trinityDE <- dexus(Trinity.trans2[,1:33], normalization ='none')
informativetranscripts <- INI(trinityDE, threshold=0.2)
plot(sort(informativetranscripts), idx = 1:30)
top30names<-(informativetranscripts@transcriptNames)[1:30]
