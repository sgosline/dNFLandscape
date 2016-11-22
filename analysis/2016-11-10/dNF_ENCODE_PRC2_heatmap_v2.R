source("../../bin/encodeSkinRNASeq.R")

##first check sanity of data from encode

#for(t in c('TPM','FPKM'))
#  for(a in c('hg19','grch38'))
#    clusterSamples(t,a)


##now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals=rna_fpkm_matrix(byIsoform=FALSE)

#count_matrix(stored=TRUE,doNorm=FALSE,minCount=2,doLogNorm=FALSE,doVoomNorm=TRUE)
skin=getGeneNamesForMatrix(getEncodeSkinMatrix(metric='FPKM',alignment='hg19',doVoomNorm = FALSE))

over=intersect(rownames(dermals),rownames(skin))
##which annotation should we do? Are they really just duplicates of one another? 

##step 1 - just combine all
comb=cbind(dermals[over,],skin[over,])

##step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr<-1:nrow(comb) #which(apply(comb,1,function(x) all(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(comb)[expr]),grep("SNO",rownames(comb)[expr]))])

##step 3, normalize
require(limma)
comb.norm=data.frame(voomWithQualityWeights(comb[expr,])$E)
comb.norm$Gene=rownames(comb.norm)

require(tidyr)
norm.df<-tidyr::gather(comb.norm,"Sample","Expression",1:(ncol(comb.norm)-1))
##prc2genes <- list("SUZ12", "EED", "EZH1", "EZH2")

EED.df <- subset(norm.df, Gene=="EED")
SUZ12.df <- subset(norm.df, Gene=="SUZ12")
EZH1.df <- subset(norm.df, Gene=="EZH1")
EZH2.df <- subset(norm.df, Gene=="EZH2")
NF1.df <- subset(norm.df, Gene=="NF1")
CDKN2A.df <- subset(norm.df, Gene=="CDKN2A")
TP53.df <- subset(norm.df, Gene=="TP53")

EED.df <- tidyr::spread(EED.df, Sample, Expression)
SUZ12.df <- tidyr::spread(SUZ12.df, Sample, Expression)
EZH1.df <- tidyr::spread(EZH1.df, Sample, Expression)
EZH2.df <- tidyr::spread(EZH2.df, Sample, Expression)
NF1.df <- tidyr::spread(NF1.df, Sample, Expression)
CDKN2A.df <- tidyr::spread(CDKN2A.df, Sample, Expression)
TP53.df <- tidyr::spread(TP53.df, Sample, Expression)

comb1.df <- dplyr::bind_rows(EED.df, SUZ12.df)
comb2.df <- dplyr::bind_rows(EZH1.df, EZH2.df)
comb3.df <- dplyr::bind_rows(NF1.df, CDKN2A.df)
comb4.df <- dplyr::bind_rows(comb1.df, comb2.df)

comb4.1 <- comb4.df[,-1]
rownames(comb4.1) <- comb4.df[,1]

comb4.2 <- as.matrix(comb4.1)

pheatmap(comb4.2, cellwidth = 10, cellheight = 15, treeheight_col= 40, treeheight_row = 40, border_color = "grey")

comb5.df <- dplyr::bind_rows(comb4.df, comb3.df)

comb5.1 <- comb5.df[,-1]
rownames(comb5.1) <- comb5.df[,1]

comb5.2 <- as.matrix(comb5.1)

pheatmap(comb5.2, cellwidth = 10, cellheight = 15, treeheight_col= 40, treeheight_row = 40, border_color = "grey")

comb6.df <- dplyr::bind_rows(comb5.df, TP53.df)

comb6.1 <- comb6.df[,-1]
rownames(comb6.1) <- comb6.df[,1]

comb6.2 <- as.matrix(comb6.1)

library(pheatmap)
pheatmap(comb6.2, cellwidth = 10, cellheight = 15, treeheight_col= 40, treeheight_row = 40, border_color = "grey")


