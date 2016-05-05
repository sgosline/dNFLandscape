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
expr<-which(apply(comb,1,function(x) all(x>0.2)))
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

##keep expresison matrix around
expr.mat<-as.matrix(comb.norm[,-ncol(comb.norm)])

##join full data frame
full.df<-join(norm.df,pat.df,by='Sample')

library(ggplot2)

##now evaluate box plots
ggplot(full.df)+geom_boxplot(mapping=aes(x=Sample,y=Expression,fill=Source))+scale_y_log10()
ggsave('allRNASeqBoxplot.png')


ggplot(subset(full.df,Gene=='NF1'))+geom_bar(mapping=aes(x=Sample,y=Expression,fill=Source),stat='identity')
ggsave('NF1_expression_sample.png')


ggplot(subset(full.df,Gene=='NF1'))+geom_boxplot(mapping=aes(x=Source,y=Expression,fill=Source))
ggsave('NF1_expression_source.png')


##now do a PC plot
pc<-prcomp(t(expr.mat),scale=T,center=T)
pcp<-ggbiplot(pc,groups=pat.df$Source,var.axes=F)
           ##now do differential expression
ggsave(pcp,file='dermalEncode_source_PCA.png')


pc<-prcomp(t(expr.mat),scale=T,center=T)
pcp<-ggbiplot(pc,groups=pat.df$Library,var.axes=F)
##now do differential expression
ggsave(pcp,file='dermalEncode_library_PCA.png')

##now just apply limma!
dm<-pat.df[,c("Source","Gender","Library")]
rownames(pat.df)<-pat.df$Sample
dm$Source<-as.character(dm$Source)
dm$Source[which(dm$Source!='dermalNF')]<-'ENCODESkin'
dm$Source<-as.factor(dm$Source)
rownames(dm)<-colnames(expr.mat)
#dm
design<-model.matrix(~Source+Gender+Library,dm)
fit <- lmFit(expr.mat, design)
fit <- eBayes(fit)
source.tab <- topTable(fit, coef='SourceENCODESkin',number=Inf,adjust.method='BY')
write.table(source.tab,file=paste('diffExValuesBetweenDermalAndAllENCODE_exprOnly.txt',sep=''))

source.tab <- topTable(fit, coef='SourceENCODESkin',number=Inf,adjust.method='BY',p.value=10e-3,lfc=1)
write.table(source.tab,quote=F,file=paste('diffExGenesBetweenDermalAndAllENCODE_10e-3_lfc1_exprOnly.txt',sep=''))
write.table(rownames(source.tab[which(source.tab$logFC<0),]),file='upRegInDermals.txt',quote=F,row.names=F,col.names=F)
write.table(rownames(source.tab[which(source.tab$logFC>0),]),file='downRegInDermals.txt',quote=F,row.names=F,col.names=F)


lib.tab <- topTable(fit, coef='LibraryRNA',number=Inf,adjust.method='BY')
write.table(lib.tab,file=paste('diffExGenesBetweenLibrary.txt_exprOnly',sep=''))

comb.tab<-topTable(fit, coef=c('SourceENCODESkin','LibraryRNA'),number=Inf,adjust.method='BY')
write.table(comb.tab,file=paste('diffExGenesBetweenDermalAndAllENCODEWithLibrary_exprOnly.txt',sep=''))

##now repeat analysis for just polyA genes? 
pa.only<-subset(dm,Library!='RNA')
pa.expr<-expr.mat[,rownames(pa.only)]
fit <- lmFit(pa.expr, model.matrix(~Source+Gender,pa.only))
fit <- eBayes(fit)
pa.tab <- topTable(fit, coef='SourceENCODESkin',number=Inf,adjust.method='BY')
write.table(pa.tab,file=paste('diffExGenesBetweenDermalAndPolyAENCODE_exprOnly.txt',sep=''))


rice.list<-c('S100B','UCHL1','ENO2','RET','GFRA2','GFRA3','ARTN','NRTN','GAP43','GDNF',
             'GFRA1','GFRA4','PSPN','NGF','BDNF','NTF4','NTF3','NTSR1','NTSR2','SORT1','NEFH','CALCA')
RiceList<-c(rep('POS',10),rep('NEG',12))
names(RiceList)<-rice.list
pheatmap(expr.mat[intersect(rice.list,rownames(expr.mat)),],annotation_col=pat.df[,c('Library','Gender','Source')],annotation_row=data.frame(RiceList),file='RiceGenesInData.png')
