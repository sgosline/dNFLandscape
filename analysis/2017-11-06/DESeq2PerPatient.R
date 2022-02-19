library(DESeq2)
## get dermal NF data
source("../../dermalNF/bin/dermalNFData.R")

this.file = "https://github.com/allaway/dNFLandscape/blob/master/analysis/2017-11-06/DESeq2PerPatient.R"

res = rna_count_matrix(stored=TRUE,doNorm=FALSE,minCount=0,doLogNorm=FALSE,doVoomNorm=FALSE)

moreids<-synTableQuery("SELECT * FROM syn5556216")@values %>% filter(!is.na(RNASeq)) %>% 
  dplyr::select(sampleIdentifier, RNASeq) %>% 
  arrange(RNASeq)

colnames(res) <- gsub("X", "", colnames(res))
colnames(res) <- as.character(gsub("\\.", "-", colnames(res)))
colnames(res) <- moreids$sampleIdentifier[colnames(res) == moreids$RNASeq]

moreids$patient <- gsub("tumor\\d+$", "", moreids$sampleIdentifier)

coldata <- moreids %>% dplyr::select(sampleIdentifier, patient) %>% 
  column_to_rownames("sampleIdentifier")

df <- as.data.frame(c())

for(i in unique(coldata$patient)){
  coldata$contrast[coldata$patient!=i] <- "ctrl"
  coldata$contrast[coldata$patient==i] <- "test"


  dds <- DESeqDataSetFromMatrix(countData = res,
                                colData = coldata,
                                design= ~ contrast)
  
  dds <- DESeq(dds)
  res.2 <- results(dds, contrast=c("contrast", "test", "ctrl"))
  res.2 <- as.data.frame(res.2)
  res.2$test <- paste0(i,"vsall")
  res.2$gene <- rownames(res.2)
  
  df <- bind_rows(df, res.2)
}

write.table(df, "OneVsAllDifferentialExpression.txt", sep ="\t")
synStore(File("OneVsAllDifferentialExpression.txt", parentId = "syn11395782"), 
         used = c("syn5051784","syn5556216"), executed = this.file)


library(ggplot2)
library(ggrepel)

ggplot(data = df, aes(x = log2FoldChange, y = -log(padj))) +
  geom_point(aes(color = padj<0.05)) +
  geom_label_repel(data = df %>% 
                     filter(padj<0.05 & abs(log2FoldChange)>10 | -log(padj) > 50), 
                   aes(x = log2FoldChange, y = -log(padj), label = gene))


df2 <- filter(df, padj < 0.05, abs(log2FoldChange) > 5)

subres<-res[rownames(res) %in% df2$gene,]
  
library(pheatmap)
pheatmap(subres)


                   