####compare TCGA expression data with dermal data

source('../../bin/TcgaExpressionData.R')
source('../../bin/TcgaMutationalData.R')

##now get dermalNF data

source("../../bin/dermalNFData.R")
rna.counts<-rna_count_matrix(doLogNorm=T,minCount=2)
#rna.fpkm<-rna_fpkm_matrix()

expr.pats<-toPatientId(colnames(alldat))
gene.symbs<-sapply(alldat[,1],function(x) unlist(strsplit(x,split='|',fixed=T))[1])
#fpkm.idx <- match(rownames(rna.fpkm),gene.symbs)
counts.idx <-match(rownames(rna.counts),gene.symbs)

disdat<-alldat[,intersect(colnames(alldat),unlist(tumsByDis))]
##can we just get all cors
full.mat<-cbind(rna.counts[!is.na(counts.idx),],disdat[counts.idx[!is.na(counts.idx)],])

pc<-prcomp(full.mat)