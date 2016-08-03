##check NF1 variants across callers
library(synapseClient)
synapseLogin()
require(dplyr)
#'get NF 1 variants from maf
getNf1RegionFromMaf<-function(maffile,synid){
  tab<-subset(read.table(gzfile(synGet(synid)@filePath),sep='\t',header=T),Hugo_Symbol=='NF1')
  sampid<-unlist(strsplit(maffile,split='.maf'))[1]
  tab$Sample<-rep(sampid,nrow(tab))
  return(tab)
}

#'vardict mafs are too big for gzfile and read.table
require(data.table)
library(R.utils)
getNf1RegionFromFullMaf<-function(maffile,synid){
   # gzf<-synGet(synid,downloadFile=F)@filePath
  #  if(!file.exists(gsub('.gz','',gzf))){
      gzf<-synGet(synid,downloadFile=T)@filePath
      gunzip(gzf)
  #  }
    tab<-subset(fread(gsub('.gz','',gzf)),Gene=='NF1')
    tab
}

##mafs are in different directories
#vardict.sid='syn6022474'
#mutect.sid='syn6186823'
#varscan.sid='syn6834373'

#
vardict.mafs<-synQuery("select Id,name from entity where parentId=='syn6022474'")
mutect.mafs<-synQuery("select Id,name from entity where parentId=='syn6186823'")
varscan.mafs<-synQuery("select Id,name from entity where parentId=='syn6834373'")
require(dplyr)
mutect.tab<-do.call('rbind',apply(mutect.mafs,1,function(x){
  getNf1RegionFromMaf(x[[1]],x[[2]])
    }))

mutect.tab$Sample<-sapply(mutect.tab$Sample,function(x) unlist(strsplit(x,split='.snp'))[1])
mutect.mods<-mutect.tab%>%filter(FILTER=='PASS')#%>%filter(IMPACT!='MODIFIER')
mutect.mods$DetectionTool=rep('Mutect',nrow(mutect.mods))

##vep failed to annotate
varscan.tab<-do.call('rbind',apply(varscan.mafs,1,function(x){
  getNf1RegionFromMaf(x[[1]],x[[2]])
}))
varscan.tab$Sample<-sapply(varscan.tab$Sample,function(x) unlist(strsplit(x,split='_VarScan'))[1])

varscan.mods<-varscan.tab%>%filter(FILTER=='PASS')#%>%filter(IMPACT!='MODIFIER')
varscan.mods$DetectionTool=rep('VarScan',nrow(varscan.mods))

vardict.tab<-do.call('rbind',apply(vardict.mafs,1,function(x){
  getNf1RegionFromFullMaf(x[[1]],x[[2]])
}))

vardict.mods<-vardict.tab%>%filter(PASS=='TRUE')%>%filter(IMPACT!='MODIFIER')
vardict.mods$DetectionTool=rep("VarDict",nrow(vardict.mods))

allmods<-rbind(varscan.mods,mutect.mods,vardict.mods)

##get dna position of mutation
allmods$DNAPos<-as.numeric(sapply(allmods$HGVSc,function(x) {
  gs<-gsub('c.','',x,fixed=T)
  gs<-gsub('*','',gs,fixed=T)
  gs2<-unlist(strsplit(gs,"N",fixed=T))[1]
  gs3<-unlist(strsplit(gs2,"-",fixed=T))[1]
  gs4<-unlist(strsplit(gs3,"+",fixed=T))[1]
  gs5<-unlist(strsplit(gs4,'_',fixed=T))[1]
  gs5}))
library(ggplot2)
#ggplot(allmods)+geom_jitter(aes(x=DNAPos,y=Sample,color=DetectionTool),width=0)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+facet_grid(. ~ IMPACT)

##now get proper patient/tumor number
tum.map<-synTableQuery('SELECT distinct Patient,TumorNumber,WGS FROM syn5556216 where WGS is not NULL')@values
allmods$PatientTumor<-sapply(allmods$Sample,function(x){
  vals<-unlist(strsplit(x,split='_'))
  pat=vals[2]
  tum<-tum.map$TumorNumber[match(vals[4],tum.map$WGS)]
  return(paste(c('Patient',pat,'Tumor',tum),collapse=' '))
})


ggplot(allmods%>%filter(IMPACT!='MODIFIER'))+geom_jitter(aes(x=DNAPos,y=PatientTumor,shape=DetectionTool,color=IMPACT))+theme(axis.text.x = element_text(angle = 90, hjust = 1))#+facet_grid(. ~ IMPACT)
ggsave("NF1SomaticMutations.png")
synStore(File("NF1SomaticMutations.png",parentId='syn6126468'))
##takeas a long time since these are full files
