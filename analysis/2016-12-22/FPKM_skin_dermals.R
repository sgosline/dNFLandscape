library(GSEABase)
library(ggplot2)
source("../../bin/encodeSkinRNASeq.R")

## now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals = rna_fpkm_matrix(byIsoform = FALSE)

# count_matrix(stored=TRUE,doNorm=FALSE,minCount=2,doLogNorm=FALSE,doVoomNorm=TRUE)
skin = getGeneNamesForMatrix(getEncodeSkinMatrix(metric = "FPKM", alignment = "hg19",
                                                 doVoomNorm = FALSE))
over = intersect(rownames(dermals), rownames(skin))

## which annotation should we do? Are they really just duplicates of one another?

## step 1 - just combine all
comb = cbind(dermals[over, ], skin[over, ])

##the PBK ids are missing from table, so need to query annotations
res<-synQuery("select patientID,tissueID,sampleID from entity where parentId=='syn5492805'")
#map<-unique(res)
#from table get generic tumor id

#tres<-synTableQuery("SELECT Patient,RnaID,TumorNumber,'RNASeq (Cufflinks)' FROM syn5556216 where RnaID is not NULL")@values
#idx<-match(res$entity.id,tres$`RNASeq (Cufflinks)`)
#dres<-res[which(!is.na(idx)),]
#tres<-tres[idx[which(!is.na(idx))],]
#full.map<-cbind(dres,tres)

this.file = 'https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2016-12-22/FPKM_skin_dermals.R'

comb = cbind("Gene" = rownames(comb), comb[,-1])
rownames(comb) <- c()
comb <- melt(comb, "Gene")
colnames(comb) <- c("Gene", "Sample", "FPKM")

comb$Type <- substr(comb$Sample, 0, 3)
comb$Type <- sub("X30", "dNF", comb$Type)
comb$Type <- sub("ENC", "skin", comb$Type)

violin<-ggplot(data=comb, aes(x=Sample, y=log2(FPKM+1))) + geom_violin(aes(color=comb$Type, fill = comb$Type))
violinall<-ggplot(data=comb, aes(x=Type, y=log2(FPKM+1))) + geom_violin(aes(color=comb$Type, fill = comb$Type))

bxplt<-ggplot(data=comb, aes(x=Sample, y=log2(FPKM+1))) + geom_boxplot(aes(color=comb$Type, fill = comb$Type))
bxpltall<-ggplot(data=comb, aes(x=Type, y=log2(FPKM+1))) + geom_boxplot(aes(color=comb$Type, fill = comb$Type))

ggsave("skin_dNF_logFPKM+1_violin.png", violin)
ggsave("skin_dNF_logFPKM+1_violinMerge.png", violinall)
ggsave("skin_dNF_logFPKM+1_boxplot.png", bxplt)
ggsave("skin_dNF_logFPKM+1_boxplotMerge.png", bxpltall)

synStore(File("skin_dNF_logFPKM+1_violin.png", parentId='syn7900368'), executed = this.file)
synStore(File("skin_dNF_logFPKM+1_violinMerge.png", parentId='syn7900368'), executed = this.file)
synStore(File("skin_dNF_logFPKM+1_boxplot.png", parentId='syn7900368'), executed = this.file)
synStore(File("skin_dNF_logFPKM+1_boxplotMerge.png", parentId='syn7900368'), executed = this.file)
