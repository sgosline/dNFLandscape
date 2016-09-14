#' Cluster cell lines by tissue type by expression
#'
#'
wd=getwd()
source("../../../RASPathwaySig/bin/cBioPortalData.R",chdir=T)
setwd(wd)
source("../../../dermalNF/bin/dermalNFData.R")
##get zscored data
load('../../../RASPathwaySig/analysis/2016-08-23/exprData.Rdata')
tcga.mat<-exprData#etDisExpressionData('',getZscores=TRUE)

#map by cell line name
#all.tiss<-unique(sapply(colnames(ccle.mat),function(x) paste(unlist(strsplit(x,split='_'))[-1],collapse='_')))
#getSamplesForDisease

tcga.dis.averages<-sapply(setdiff(tcga.cancer.types,c('lgggbm','nsclc')),function(x) {
  samps<-getSamplesForDisease(x)
  samps<-sapply(samps,function(y) gsub('-','.',y,fixed=T))
  cols<-match(samps,colnames(tcga.mat))
  cols<-cols[!is.na(cols)]
#  cols=grep(x,colnames(ccle.mat))
  if(length(cols)>1)
    return(rowMeans(tcga.mat[,cols],na.rm=T))
  else
    return(tcga.mat[,cols])})

##get pnf cell line data
require(synapseClient)
synapseLogin()

#get the TPM values

zscore<-function(x){
  x<-unlist(x)
  (x-mean(x,na.rm=T))/sd(x)
}


dermalSamps<-rna_fpkm_matrix()
#read.table(synGet('syn5580347')@filePath,sep='\t',header=T)

#z-score by gene
normDermData<-apply(dermalSamps,2,zscore)

phenoData<-fpkm_annotations()
phenoData<-unique(phenoData[,c(1,2,5)])
pats<-paste('patient',apply(phenoData[match(colnames(normDermData),phenoData$sample),c(1,3)],1,paste,collapse=' tumor'))

colnames(normDermData)<-pats


##get column names

##get other nf data
comm.genes<-intersect(rownames(tcga.dis.averages),rownames(normDermData))
library(pheatmap)
cmat<-cbind(normDermData[comm.genes,],tcga.dis.averages[comm.genes,])

##let's plot!
png('tcga_dermal_dendrogram.png',width=800)
plot(hclust(dist(t(cmat))),main='TCGA Samples with dermal NF data')
dev.off()

#synStore(File('tcga_dermal_dendrogram.png',parentId='syn5594111'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2016-08-12/clusterCellLinesByExpression.R')),used=list(list(entity='syn5580347')))

png('tcga_dermal_cor_dendrogram.png',width=800)
plot(hclust(as.dist(1-cor(cmat,use='pairwise.complete.obs'))),main='TCGA Samples with Dermal Data')
dev.off()

library(ggbiplot)
rcol<-setdiff(colnames(cmat),c('pcpg','lihc'))
pc<-prcomp(t(cmat),center=T,scale=T)
#synStore(File('ccle_pnf_dendrogram.pdf',parentId='syn5594111'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2016-08-12/clusterCellLinesByExpression.R')),used=list(list(entity='syn5580347')))

ggbiplot(pc,var.axes=F,labels=colnames(cmat),choices=1:2,center=T,scale=T)
ggsave('tcga_dermalNF_samps_pca12.png')
ggbiplot(pc,var.axes=F,labels=colnames(cmat),choices=2:3,center=T,scale=T)
ggsave('tcga_dermalNF_samps_pca23.png')

for(f in list.files('.','.png'))
  synStore(File(f,parentId='syn5821631'),
           executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-09-06/clusterCellLinesByExpression.R')))

#synStore(File('ccle_pnf_cor_dendrogram.pdf',parentId='syn5594111'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/pnfCellLines/master/analysis/2016-08-12/clusterCellLinesByExpression.R')),used=list(list(entity='syn5580347')))
