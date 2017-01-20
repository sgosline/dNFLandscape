library(ggplot2)
library(data.table)
library(reshape2)
library(synapseClient)

rnaid<-synTableQuery('SELECT RNASeq, sampleIdentifier FROM syn5556216')@values
rnaid<-dplyr::filter(rnaid, !is.na(RNASeq))
rnaid$RNASeq<-lapply(rnaid$RNASeq, function(x) synGet(id=x, downloadFile = FALSE)@fileHandle$fileName)
rnaid$RNASeq<-gsub('accepted_hits_', "", rnaid$RNASeq)
rnaid$RNASeq<-gsub('_featureCounts.txt', "", rnaid$RNASeq)

TrinityTMM.genes <- synGet('syn8034122')@filePath
Trinity.genes <- fread(file = TrinityTMM.genes, sep = '\t', header = TRUE)
colnames(Trinity.genes)<-rnaid$sampleIdentifier[match(colnames(Trinity.genes), rnaid$RNASeq)]
Trinity.genes <- melt(Trinity.genes)
colnames(Trinity.genes) <- c("Gene", "Sample", "TMM")

TrinityTMM.transcripts <- synGet('syn8034127')@filePath
Trinity.transcripts <- fread(file = TrinityTMM.transcripts, sep = '\t', header = TRUE)
colnames(Trinity.transcripts)<-rnaid$sampleIdentifier[match(colnames(Trinity.transcripts), rnaid$RNASeq)]
Trinity.transcripts <- melt(Trinity.transcripts)
colnames(Trinity.transcripts) <- c("Transcript", "Sample", "TMM")

this.file = 'https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2017-01-05/Trinity_TMM.R'

violin.genes<-ggplot(data=Trinity.genes, aes(x=Sample, y=log2(TMM+1))) + geom_violin(aes(color=Sample, fill = Sample))
violin.transcripts<-ggplot(data=Trinity.transcripts, aes(x=Sample, y=log2(TMM+1))) + geom_violin(aes(color=Sample, fill = Sample))
bxplt.genes<-ggplot(data=Trinity.genes, aes(x=Sample, y=log2(TMM+1))) + geom_boxplot(aes(color=Sample, fill = Sample))
bxplt.transcripts<-ggplot(data=Trinity.transcripts, aes(x=Sample, y=log2(TMM+1))) + geom_boxplot(aes(color=Sample, fill = Sample))

ggsave("Trinity_genes_log2TMM+1_violin_SECONDBUILD.png", violin.genes)
ggsave("Trinity_transcripts_log2TMM+1_violin_SECONDBUILD.png", violin.transcripts)
ggsave("Trinity_genes_log2TMM+1_boxplot_SECONDBUILD.png", bxplt.genes)
ggsave("Trinity_transcripts_log2TMM+1_boxplot_SECONDBUILD.png", bxplt.transcripts)

synStore(File("Trinity_genes_log2TMM+1_violin_SECONDBUILD.png", parentId='syn8027359'), executed = this.file, used = 'syn7989909')
synStore(File("Trinity_transcripts_log2TMM+1_violin_SECONDBUILD.png", parentId='syn8027359'), executed = this.file, used = 'syn7989913')
synStore(File("Trinity_genes_log2TMM+1_boxplot_SECONDBUILD.png", parentId='syn8027359'), executed = this.file, used = 'syn7989909')
synStore(File("Trinity_transcripts_log2TMM+1_boxplot_SECONDBUILD.png", parentId='syn8027359'), executed = this.file, used = 'syn7989913')


