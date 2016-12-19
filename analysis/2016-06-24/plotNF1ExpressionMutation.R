source('../../../dermalNF/bin/dermalNFData.R')
source("../../../dermalNF/bin/WGSData_VarDict.R")


nf1.analysis.dir='syn6175555'

rna.seq<-rna_fpkm_matrix(TRUE)
nf1.mut<-getMutationStatsForGene(all.gene.muts,gene='NF1')

nf1.expr<-rna.seq[grep('^NF1.NM',rownames(rna.seq)),]

trans<-sapply(rownames(nf1.expr),function(x) unlist(strsplit(x,split='.',fixed=T))[2])

idx<-match(sapply(as.character(nf1.mut$Transcript),function(x) unlist(strsplit(x,split='.',fixed=T))[1]),trans)


##now we have to match the patient tumors!!!