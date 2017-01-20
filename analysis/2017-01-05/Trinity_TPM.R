library(ggplot2)
library(data.table)
library(reshape2)
library(synapseClient)

rnaid<-synTableQuery('SELECT RNASeq, sampleIdentifier FROM syn5556216')@values
rnaid<-dplyr::filter(rnaid, !is.na(RNASeq))
rnaid$RNASeq<-lapply(rnaid$RNASeq, function(x) synGet(id=x, downloadFile = FALSE)@fileHandle$fileName)
rnaid$RNASeq<-gsub('accepted_hits_', "", rnaid$RNASeq)
rnaid$RNASeq<-gsub('_featureCounts.txt', "", rnaid$RNASeq)

TrinityTPM.genes <- synGet('syn7989908')@filePath
Trinity.genes <- fread(file = TrinityTPM.genes, sep = '\t', header = TRUE)
colnames(Trinity.genes)<-rnaid$sampleIdentifier[match(colnames(Trinity.genes), rnaid$RNASeq)]
Trinity.genes <- melt(Trinity.genes)
colnames(Trinity.genes) <- c("Gene", "Sample", "TPM")

TrinityTPM.transcripts <- synGet('syn7989912')@filePath
Trinity.transcripts <- fread(file = TrinityTPM.transcripts, sep = '\t', header = TRUE)
colnames(Trinity.transcripts)<-rnaid$sampleIdentifier[match(colnames(Trinity.transcripts), rnaid$RNASeq)]
Trinity.transcripts <- melt(Trinity.transcripts)
colnames(Trinity.transcripts) <- c("Transcript", "Sample", "TPM")

this.file = 'https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2017-01-05/Trinity_TMM.R'

violin.genes<-ggplot(data=Trinity.genes, aes(x=Sample, y=log2(TPM+1))) + geom_violin(aes(color=Sample, fill = Sample))
violin.transcripts<-ggplot(data=Trinity.transcripts, aes(x=Sample, y=log2(TPM+1))) + geom_violin(aes(color=Sample, fill = Sample))
bxplt.genes<-ggplot(data=Trinity.genes, aes(x=Sample, y=log2(TPM+1))) + geom_boxplot(aes(color=Sample, fill = Sample))
bxplt.transcripts<-ggplot(data=Trinity.transcripts, aes(x=Sample, y=log2(TPM+1))) + geom_boxplot(aes(color=Sample, fill = Sample))

ggsave("Trinity_genes_log2TPM+1_violin.png", violin.genes)
ggsave("Trinity_transcripts_log2TPM+1_violin.png", violin.transcripts)
ggsave("Trinity_genes_log2TPM+1_boxplot.png", bxplt.genes)
ggsave("Trinity_transcripts_log2TPM+1_boxplot.png", bxplt.transcripts)

synStore(File("Trinity_genes_log2TPM+1_violin.png", parentId='syn7986689'), executed = this.file, used = 'syn7989908')
synStore(File("Trinity_transcripts_log2TPM+1_violin.png", parentId='syn7986689'), executed = this.file, used = 'syn7989912')
synStore(File("Trinity_genes_log2TPM+1_boxplot.png", parentId='syn7986689'), executed = this.file, used = 'syn7989908')
synStore(File("Trinity_transcripts_log2TPM+1_boxplot.png", parentId='syn7986689'), executed = this.file, used = 'syn7989912')

Trinity.genes.m <- as.matrix(fread(file = TrinityTPM.genes, sep = '\t', header = TRUE))
library(pheatmap)
pheatmap(Trinity.genes.m)
