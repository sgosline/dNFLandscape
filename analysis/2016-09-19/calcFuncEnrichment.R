##
source('../../dermalNF/bin/dermalNFData.R')

##first let's figure out which WGCNA run to use.
quant=c("cuffLinks",'featureCounts')
filter=c("Conn",'Sd')
sign=c('signed','unsigned')
pref='TOM'
suffix='top10GOTermsPermodule.csv'

allfiles=synQuery('select name,id from entity where parentId=="syn5669860"')

muts<-read.table(synGet('syn6097853')@filePath,sep='\t',header=T)


library(dplyr)
res<-c()
for(q in quant){
  for(f in filter){
    for(s in sign){
      fname=paste(pref,q,'_','filterBy',f,'_',s,suffix,sep='')
      eid=allfiles$entity.id[match(fname,allfiles$entity.name)]
      if(is.na(eid)){
        fname=paste(pref,q,'_','filteredBy',f,'_',s,suffix,sep='')
        eid=allfiles$entity.id[match(fname,allfiles$entity.name)]
      }
      tab<-read.table(synGet(eid)@filePath,sep=',',header=T,as.is=T)
      sigs<-subset(tab,BonferoniP<0.1)
      num.modules<-length(unique(sigs$module))
      num.terms<-length(unique(sigs$termName))
      tres<-sigs%>%group_by(module)%>%summarize(TPM=n_distinct(termName))
      terms.perm.module<-mean(tres$TPM)
      res=rbind(res,list(Quant=q,Filter=f,Sign=s,SigTerms=nrow(sigs),NumModules=num.modules,NumTerms=num.terms,
                     TermsPerModule=terms.perm.module))
      
    
    }
  }
}

View(res)

##now figure out samples in each cluster. 


