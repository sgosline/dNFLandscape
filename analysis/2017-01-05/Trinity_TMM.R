library(ggplot2)
library(data.table)
library(reshape2)
library(synapseClient)

TrinityTMM.genes <- synGet('syn7989909')@filePath
Trinity.genes <- fread(file = TrinityTMM.genes, sep = '\t', header = TRUE)
Trinity.genes <- melt(Trinity.genes)
colnames(Trinity.genes) <- c("Gene", "Sample", "TMM")

TrinityTMM.transcripts <- synGet('syn7989913')@filePath
Trinity.transcripts <- fread(file = TrinityTMM.transcripts, sep = '\t', header = TRUE)
Trinity.transcripts <- melt(Trinity.transcripts)
colnames(Trinity.transcripts) <- c("Transcript", "Sample", "TMM")

this.file = 'https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2017-01-05/Trinity_TMM.R'

violin.genes<-ggplot(data=Trinity.genes, aes(x=Sample, y=log2(TMM+1))) + geom_violin(aes(color=Sample, fill = Sample))
violin.transcripts<-ggplot(data=Trinity.transcripts, aes(x=Sample, y=log2(TMM+1))) + geom_violin(aes(color=Sample, fill = Sample))
bxplt.genes<-ggplot(data=Trinity.genes, aes(x=Sample, y=log2(TMM+1))) + geom_boxplot(aes(color=Sample, fill = Sample))
bxplt.transcripts<-ggplot(data=Trinity.transcripts, aes(x=Sample, y=log2(TMM+1))) + geom_boxplot(aes(color=Sample, fill = Sample))

ggsave("Trinity_genes_log2TMM+1_violin.png", violin.genes)
ggsave("Trinity_transcripts_log2TMM+1_violin.png", violin.transcripts)
ggsave("Trinity_genes_log2TMM+1_boxplot.png", bxplt.genes)
ggsave("Trinity_transcripts_log2TMM+1_boxplot.png", bxplt.transcripts)

synStore(File("Trinity_genes_log2TMM+1_violin.png", parentId='syn7986689'), executed = this.file, used = 'syn7989909')
synStore(File("Trinity_transcripts_log2TMM+1_violin.png", parentId='syn7986689'), executed = this.file, used = 'syn7989913')
synStore(File("Trinity_genes_log2TMM+1_boxplot.png", parentId='syn7986689'), executed = this.file, used = 'syn7989909')
synStore(File("Trinity_transcripts_log2TMM+1_boxplot.png", parentId='syn7986689'), executed = this.file, used = 'syn7989913')


