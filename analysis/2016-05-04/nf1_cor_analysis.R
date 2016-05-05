##nf1 correlation studies

source("../../bin/encodeSkinRNASeq.R")

##first check sanity of data from encode

#for(t in c('TPM','FPKM'))
#  for(a in c('hg19','grch38'))
#    clusterSamples(t,a)


##now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals=rna_fpkm_matrix(doVoomNorm=FALSE)

#count_matrix(stored=TRUE,doNorm=FALSE,minCount=2,doLogNorm=FALSE,doVoomNorm=TRUE)
skin=getGeneNamesForMatrix(getEncodeSkinMatrix(metric='FPKM',alignment='hg19',doVoomNorm = FALSE))

over=intersect(rownames(dermals),rownames(skin))
##which annotation should we do? Are they really just duplicates of one another? 

##step 1 - just combine all
comb=cbind(dermals[over,],skin[over,])

##step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr<-which(apply(comb,1,function(x) any(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(comb)[expr]),grep("SNO",rownames(comb)[expr]))])



##step 3, normalize
require(limma)
comb.norm=data.frame(voomWithQualityWeights(comb[expr,])$E)
comb.norm$Gene=rownames(comb.norm)

require(tidyr)
norm.df<-tidyr::gather(comb.norm,"Sample","Expression",1:(ncol(comb.norm)-1))

samps<-getSampleNamesForMatrix(skin)
dermal.vars<-rna_patient_variables()
dinds<-sapply(colnames(dermals),function(x) match(gsub("X","",gsub('.','-',x,fixed=T)),dermal.vars$sampleID))
dermal.vars<-dermal.vars[dinds,]

##full patient annotation
pat.df<-data.frame(Source=c(rep('dermalNF',ncol(dermals)),as.character(samps$Sample)),
                   Gender=c(tolower(as.character(dermal.vars$Gender)),as.character(samps$Sex)),
                   Age=c(as.character(dermal.vars$Age),sapply(as.character(samps$Age),function(x) gsub(' years','',x))),
                   Sample=c(colnames(dermals),c(colnames(skin))),
                   Library=c(rep('polyadenylated mRNA',ncol(dermals)),as.character(samps$Library)))

rownames(pat.df)<-pat.df$Sample
##keep expresison matrix around
expr.mat<-as.matrix(comb.norm[,-ncol(comb.norm)])

##join full data frame
full.df<-join(norm.df,pat.df,by='Sample')

fisher.z<-function(x) {(0.5)*log((1+x)/(1-x))}

nf1.cors<-apply(expr.mat,1,function(x) cor(x,expr.mat['NF1',]))

nf1.cors<-sort(nf1.cors,dec=T)

nf1.z<-sapply(nf1.cors,fisher.z)

nf1.pvals<-pnorm(nf1.z,mean=mean(nf1.z[-1]),sd=sd(nf1.z[-1]),lower.tail=F)
write.table(names(which(nf1.pvals<0.05)),file='nf1corp05.txt',row.names=F,col.names=F,quote=F)
write.table(names(which(nf1.pvals<0.01)),file='nf1corp01.txt',row.names=F,col.names=F,quote=F)
pheatmap(expr.mat[names(which(nf1.pvals<0.005)),],annotation_col=pat.df[,c('Source','Library','Gender')],cellwidth=10,cellheight=10,file='nf1CorrelatedGenesAcrossSamples.png')

nf1Correlation<-function(full.df,gene='BRAF'){
  res<-tidyr::spread(subset(full.df,Gene%in%c('NF1',gene)),Gene,Expression)
  cv<-cor(res[,'NF1'],res[,gene])
  ggplot(res)+geom_point(aes_string(x='NF1',y=gene,col='Source'))+ggtitle(paste("R=",format(cv,3)))
}

