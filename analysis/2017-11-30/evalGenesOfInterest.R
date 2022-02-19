##here we collect a list of genes and evaluate their expression in the proteomics and RNA samples

genelist<-c('RET','GFRA1','GDNF','NRTN','ARTN','PSPN','NGFR',
            'GAP43','NTRK1','NTRK2','NTRK3','NGF','BDNF','NTF3','NTF4',
            'UCHL1','S100B','ENO2','GFRA2','GFRA3','GFRA4','TGFB1','TGFB2')

#added 4/27/18
genelist <- c(genelist,'L1CAM', 'SOX10', 'GFAP', 'NCAM1', 'NRG1')

#added 5/8
genelist <- c(genelist,'TRPV1','TRPA1')

##added 8/14
genelist <-c(genelist,'AIF1','SOX2','PRSS56')

source("../../../dermalNF/bin/dermalNFData.R")

this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2017-11-30/evalGenesOfInterest.R'
#get rna data
rna.annot<-rna_annotations()
rna<-rna_count_matrix(doLogNorm=TRUE,minCount=5)

allvars=apply(rna,1,function(x) var(x,na.rm=T)^2/mean(x,na.rm=T))
pats<-rna.annot$individualID
names(pats)<-rna.annot$id

library(pheatmap)
#pheatmap(rna[order(allvars,decreasing=T)[1:100],],annotation_col=data.frame(Patient=pats),
#         cellwidth=10,cellheight=10,
#          clustering_distance_cols = 'correlation',clustering_distance_rows = 'correlation',
#         file='RNA_100_mostVariable_min5counts.png')

colnames(rna)<-paste(pats,rna.annot$specimenID)
#plot(hclust(dist(t(rna)),method='ward.D2'))


#now get all mrnAs, even the poorly expressed
rna<-rna_count_matrix(doNorm=TRUE,minCount=0)

red.set<-rna[match(genelist,rownames(rna)),]
pheatmap(log10(0.00001+red.set),annotation_col=data.frame(Patient=pats),cellwidth=10,cellheight=10,
         clustering_distance_cols = 'correlation',file='signalingProteinGeneExpress.pdf')
synStore(File('signalingProteinGeneExpress.pdf',parentId='syn6242409'),executed=list(this.script))

##now let's try ranked plost
pats.only<-rna[,names(pats)[which(pats%in%c("patient2",'patient11'))]]

ranked<-data.frame(apply(rna,2,function(x) rank(x)/length(x)))#function(x) length(x)-rank(x)))
ranked$Gene<-rownames(ranked)
require(tidyverse)
res<-tidyr::gather(ranked,key='id',value='Rank',-Gene)%>%
  left_join(select(rna.annot,id,specimenID,individualID),by="id")
newres<-res%>%mutate(inList=Gene%in%genelist)


require(ggplot2)

newres$Sample <- sapply(newres$specimenID,function(x) gsub("patient","Patient ",gsub("tumor"," Tumor ",x)))
newres$Patient <-sapply(newres$Sample,function(x) gsub(" Tumor [0-9]*","",x))

genesOfInterest<-subset(newres,inList==TRUE)

gene.means<-genesOfInterest%>%group_by(Gene)%>%summarize(meanRank=mean(Rank))

genesOfInterest$Gene<-factor(genesOfInterest$Gene,levels=gene.means$Gene[order(gene.means$meanRank)])

#create quantiles

p2<-ggplot(genesOfInterest)+geom_point(aes(x=Gene,y=Rank,col=Patient))+theme(axis.text.x=element_text(angle = -90, hjust = 0))+ggtitle('Rank of Selected genes in cNF patient samples')
ggsave('AllGeneExpressionRank.pdf')
synStore(File('AllGeneExpressionRank.pdf',parentId='syn6242409'),executed=list(this.script))
genesOfInterest<-subset(genesOfInterest,Patient%in%c('Patient 2','Patient 11'))
gene.means<-genesOfInterest%>%group_by(Gene)%>%summarize(meanRank=mean(Rank))
genesOfInterest$Gene<-factor(genesOfInterest$Gene,levels=gene.means$Gene[order(gene.means$meanRank)])


p2<-ggplot(genesOfInterest)+geom_point(aes(x=Gene,y=Rank,col=Sample,shape=Patient))+theme(axis.text.x=element_text(angle = -90, hjust = 0))+ggtitle('Rank of Selected genes in cNF patient samples')
ggsave('AllGeneExpressionRank2pats.pdf')

synStore(File('AllGeneExpressionRank2pats.pdf',parentId='syn6242409'),executed=list(this.script))

##plot all patients
mat<-spread(res,'Gene',value='Rank')%>%select(-c(id,individualID))
rownames(mat)<-mat$specimenID
mat<-select(mat,-specimenID)
gene.cors<-apply(mat,2,function(x) cor(x,mat[,'S100B']))
gene.corts<-apply(mat,2,function(x) cor.test(x,mat[,'S100B'])$p.value)