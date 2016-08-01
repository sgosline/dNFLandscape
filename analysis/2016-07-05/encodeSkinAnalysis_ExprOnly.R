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
expr<-1:nrow(comb) #which(apply(comb,1,function(x) all(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(comb)[expr]),grep("SNO",rownames(comb)[expr]))])



##step 3, normalize
#require(limma)
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
#require
#full.df<-join(norm.df,pat.df,by='Sample')

library(ggplot2)


##now just apply limma!
# dm<-pat.df[,c("Source","Gender","Library")]
 rownames(pat.df)<-pat.df$Sample
# dm$Source<-as.character(dm$Source)
# dm$Source[which(dm$Source!='dermalNF')]<-'ENCODESkin'
# dm$Source<-as.factor(dm$Source)
# rownames(dm)<-colnames(expr.mat)
# #dm
# design<-model.matrix(~Source+Gender+Library,dm)
# fit <- lmFit(expr.mat, design)
# fit <- eBayes(fit)
# #source.tab <- topTable(fit, coef='SourceENCODESkin',number=Inf,adjust.method='BY')
# 
# ##now repeat analysis for just polyA genes? 
# pa.only<-subset(dm,Library!='RNA')
# pa.expr<-expr.mat[,rownames(pa.only)]
# fit <- lmFit(pa.expr, model.matrix(~Source+Gender,pa.only))
# fit <- eBayes(fit)
# pa.tab <- topTable(fit, coef='SourceENCODESkin',number=Inf,adjust.method='BY')
# write.table(pa.tab,file=paste('diffExGenesBetweenDermalAndPolyAENCODE_exprOnly.txt',sep=''))
require(dplyr)
marco.list<-read.csv('NF1_NF2_TargetsFromLitv1.csv',header=T)

synStore(File('NF1_NF2_TargetsFromLitv1.csv',parentId='syn6242409'))


non.blank<-subset(marco.list,GeneName!="")[,1:4]

with.nf<-marco.list %>% mutate(NF1=Disease=='NF1',NF2=Disease=="NF2") %>% select(GeneName,Symptom,NF1,NF2)

###should automate this, but too much work at this point
all.types<-with.nf%>% mutate(MPNST=Symptom=="MPNST",Pain=Symptom=="PAIN",
                             OPG=Symptom=="OPG",LGG=Symptom=="LGG",Bone=Symptom=="BONE",
                             Cognition=Symptom=="COGNITION",Meningioma=Symptom=="MENINGIOMA",
                             PylocticAstrocytoma=Symptom=="PYLOCTIC ASTROCYTOMA")%>% select(-matches('Symptom'))

df<-unique(all.types)
new.df<-t(sapply(unique(df$GeneName),function(x) apply(subset(df,GeneName==x)[,-1],2,function(y) as.character(any(y)))))
colnames(new.df)<-colnames(df)[-1]
rownames(new.df)<-unique(df$GeneName)
new.df<-data.frame(new.df)
library(pheatmap)
library(ggplot2)
#marcoList<-c(rep('POS',10),rep('NEG',12))
#names(marcoList)<-.list
pheatmap(expr.mat[intersect(rownames(new.df),rownames(expr.mat)),],cellwidth = 10,cellheight=10,
         annotation_col=select(pat.df,Source),
         annotation_row=data.frame(new.df),filename='NF12GenesInNFAndSkinData.png',height=12)

nf1.genes<-subset(new.df,NF1==TRUE) %>% select(MPNST,Pain,OPG,LGG,Bone,Cognition,PylocticAstrocytoma)
pheatmap(expr.mat[intersect(rownames(nf1.genes),rownames(expr.mat)),],cellwidth = 10,cellheight=10,
         annotation_col=select(pat.df,Source),
         annotation_row=nf1.genes,file='NF1_only_GenesInNFAndSkinData.png',height=10)

nf2.genes<-subset(new.df,NF2==TRUE) %>% select(Meningioma)
pheatmap(expr.mat[intersect(rownames(nf2.genes),rownames(expr.mat)),],cellwidth = 10,cellheight=10,
         annotation_col=select(pat.df,Source),
         annotation_row=nf2.genes,file='NF2_only_GenesInNFAndSkinData.png',height=10)

##now let's do the same for the dermal only data
norm.dermals=rna_fpkm_matrix(doVoomNorm=TRUE)

pheatmap(norm.dermals[intersect(rownames(new.df),rownames(norm.dermals)),],cellwidth = 10,cellheight=10,
        # annotation_col=select(pat.df,Source),
         annotation_row=data.frame(new.df),filename='NF12GenesInNFOnlyData.png',height=12)


nf1.genes<-subset(new.df,NF1==TRUE) %>% select(MPNST,Pain,OPG,LGG,Bone,Cognition,PylocticAstrocytoma)
pheatmap(norm.dermals[intersect(rownames(nf1.genes),rownames(norm.dermals)),],cellwidth = 10,cellheight=10,
       #  annotation_col=select(pat.df,Source),
         annotation_row=nf1.genes,file='NF1_only_GenesInNFOnlyData.png',height=10)

nf2.genes<-subset(new.df,NF2==TRUE) %>% select(Meningioma)
pheatmap(norm.dermals[intersect(rownames(nf2.genes),rownames(norm.dermals)),],cellwidth = 10,cellheight=10,
       #  annotation_col=select(pat.df,Source),
         annotation_row=nf2.genes,file='NF2_only_GenesInNFOnlyData.png',height=10)

for(fi in list.files('./')[grep('png',list.files('./'))]){
  synStore(File(fi,parentId='syn6242409'),used=list(list(entity='syn6242733')),
           executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-07-05/encodeSkinAnalysis_ExprOnly.R')))
}
