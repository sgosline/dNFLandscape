source("../../bin/encodeSkinRNASeq.R")

##now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals=rna_fpkm_matrix(byIsoform=FALSE)

##step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr<-1:nrow(dermals) #which(apply(comb,1,function(x) all(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(dermals)[expr]),grep("SNO",rownames(dermals)[expr]))])

##step 3, normalize
require(limma)
dermals.norm=data.frame(voomWithQualityWeights(dermals[expr,])$E)
dermals.norm$Gene=rownames(dermals.norm)

require(dplyr)
dermals.norm <- select(dermals.norm, -Gene)

dermals.m <- as.matrix(dermals.norm)

