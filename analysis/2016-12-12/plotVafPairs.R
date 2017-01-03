source("../../bin/wgsAnalysis.R")
library(ggplot2)
library(tidyr)
##now take allele frequency for each pair for each patient

samples<-as.character(unique(expr.gene.muts1$Sample))

patients<-sapply(samples,function(x) unlist(strsplit(x,split='_'))[2])


plotVafByPairs<-function(sampMat,patient,mutDf){
    somDf<-subset(mutDf,Status%in%c('StrongSomatic','LikelySomatic'))#,'LikelyLOH'))
    somDf<-subset(somDf,!Chr%in%c('chrX','chrY'))
    pdf(paste('patient',patient,'varAlleleFreqBySamp.pdf',sep='_'))
    apply(sampMat,2,function(x){
        pairs<-subset(somDf,Sample%in%x)
        if(length(unique(pairs$Sample))==1)
          return(NULL)
        #res<-pairs%>%select(Sample,Gene,AlleleFreq)
        res<-pairs%>%unite(Pos,Chr,Start,sep='_')%>%select(Sample,Gene,Pos,AlleleFreq,VD,RD)
       # mafs<-unique(res)%>%spread(AlleleFreq,Sample,fill=0.0)
        amat<-acast(res,Pos~Sample,fill=0.0,value.var='AlleleFreq',fun.aggregate=mean)
        newdf<-data.frame(amat,Gene=res$Gene[match(rownames(amat),res$Pos)])
        gz<-which(apply(newdf[,1:2],1,function(x) all(x>0)))
        p<-ggplot(newdf,mapping=aes_string(x=x[1],y=x[2]))+geom_point()+geom_text(data=newdf[gz,],mapping=aes_string(x=x[1],y=x[2],label="Gene"))
        print(p)
          })
    dev.off()
}
sapply(unique(patients),function(pat){
  samps<-samples[which(patients==pat)]
  pars<-combn(samps,2)
  plotVafByPairs(pars,patient=pat,expr.gene.muts1)
  })
  #print(dim(pars))
for(file in list.files('.'))
  if(length(grep('pdf',file)>0))
    synStore(File(file,parentId=''),executed=list(list(executed=this.script)),used=list(list(entity='')))


