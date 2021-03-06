source("../../bin/wgsAnalysis.R")
library(ggplot2)
library(tidyr)
library(data.table)
require(parallel)
##now take allele frequency for each pair for each patient
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-19/plotVafPairsWithAll.R'

samples<-as.character(unique(expr.gene.muts1$Sample))

patients<-sapply(samples,function(x) unlist(strsplit(x,split='_'))[2])

pfiles<-synapseQuery("select * from entity where parentId=='syn5605256'")

getFileByPatient<-function(x){
  pfiles[grep(paste('pval0.05patient',x,'impact',sep='_'),pfiles$entity.name),]$entity.id
  
}
plotVafByPairs<-function(sampMat,patient,vdFilter=0){
  mutDf<-as.data.table(fread(synGet(getFileByPatient(patient))@filePath,sep='\t'))
  
    somDf<-subset(mutDf,Status%in%c('StrongSomatic','LikelySomatic'))#,'LikelyLOH'))
    somDf<-subset(somDf,!Chr%in%c('chrX','chrY'))
    somDf<-subset(somDf,VD>vdFilter)
    
    pdf(paste('patient',patient,'varAlleleFreqWithDepth_gt',vdFilter,'BySamp.pdf',sep='_'))
    sapply(unique(c(sampMat)),function(x){
      plot(density(subset(somDf,Sample==x)$AlleleFreq),main=paste("Allele Frequency for",x))
    })
    apply(sampMat,2,function(x){
        pairs<-subset(somDf,Sample%in%x)
        if(length(unique(pairs$Sample))==1)
          return(NULL)
        #res<-pairs%>%select(Sample,Gene,AlleleFreq)
        res<-pairs%>%unite(Pos,Chr,Start,sep='_')%>%select(Sample,Gene,Pos,AlleleFreq,VD,RD,Status,Effect)
               # mafs<-unique(res)%>%spread(AlleleFreq,Sample,fill=0.0)
        amat<-acast(res,Pos~Sample,fill=0.0,value.var='AlleleFreq',fun.aggregate=mean)
        newdf<-data.frame(amat,Gene=res$Gene[match(rownames(amat),res$Pos)],Status=res$Status[match(rownames(amat),res$Pos)],Effect=res$Effect[match(rownames(amat),res$Pos)])
        gz<-intersect(which(apply(newdf[,1:2],1,function(x) all(x>0))),which(newdf$Effect!='MODIFIER'))
        p<-ggplot(newdf,mapping=aes_string(x=x[1],y=x[2]))+geom_point(aes(color=Status))+geom_text(data=newdf[gz,],mapping=aes_string(x=x[1],y=x[2],label="Gene"))
        print(p)
          })
    dev.off()
}
mclapply(unique(patients),function(pat){
  samps<-samples[which(patients==pat)]
  pars<-combn(samps,2)
 # plotVafByPairs(pars,patient=pat,vdFilter=0)
#  plotVafByPairs(pars,patient=pat,vdFilter=10)
#  plotVafByPairs(pars,patient=pat,vdFilter=20)
  plotVafByPairs(pars,patient=pat,vdFilter=50)
  },mc.cores=3)
  #print(dim(pars))
#for(file in list.files('.'))
#  if(length(grep('pdf',file)>0))
#    synStore(File(file,parentId=''),executed=list(list(executed=this.script)),used=list(list(entity='')))


