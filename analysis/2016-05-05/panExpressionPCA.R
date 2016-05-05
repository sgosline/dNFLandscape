####compare TCGA expression data with dermal data

source('../../bin/TcgaExpressionData.R')
source('../../bin/TcgaMutationalData.R')

##now get dermalNF data

source("../../dermalNF/bin/dermalNFData.R")
rna.counts<-rna_count_matrix(doLogNorm=F,minCount=0)
#rna.fpkm<-rna_fpkm_matrix()

expr.pats<-toPatientId(colnames(alldat))
gene.symbs<-sapply(alldat[,1],function(x) unlist(strsplit(x,split='|',fixed=T))[1])
#fpkm.idx <- match(rownames(rna.fpkm),gene.symbs)
counts.idx <-match(rownames(rna.counts),gene.symbs)

disdat<-alldat[,intersect(colnames(alldat),unlist(tumsByDis))]
##can we just get all cors
full.mat<-cbind(rna.counts[!is.na(counts.idx),],disdat[counts.idx[!is.na(counts.idx)],])

require(tidyr)
full.df<-data.frame(full.mat)
full.df$Gene<-rownames(full.mat)
mat.df<-tidyr::gather(full.df,"SampleID","ExpressionValue",1:(ncol(full.df)-1))

disease<-sapply(mat.df$SampleID,function(x) if(length(grep("syn",x))>0) return("dermalNF")
                else return names(tumsByDis)[grep(x,tumsByDis)])


require(ggplot2)
ggplot(mat.df)+geom_boxplot(aes(y=ExpressionValue,x=SampleID))

require(limma)
##now try to do some normalization
res=voomWithQualityWeights(full.mat,design='')

##now reformat and plot to see if it worked

pc<-prcomp(full.mat)