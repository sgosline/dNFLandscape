###goal is to download dermal somatic mutations and run networks

require(synapseClient)

synapseLogin()
require(dplyr)

##get all cancer mutations
#allCancerMuts=synGet('syn5611520')

##maybe we don't want filtered, want to include common variants!
#allMuts <- synGet("syn5839666")
#allMuts <-synGet("syn5713423")
##get all mutations
source("../../bin/wgsAnalysis.R")

##get expressed genes!
#source("../../bin/dermalNFData.R")
#rna.ex=rna_count_matrix(minCount=3)
#expressed.genes=rownames(rna.ex)
#print(paste('Found',length(expressed.genes),'and reducing mutations accordingly'))


getCountsFromTable<-function(tab,mutationType=c('germline','somatic'),effect="HIGH_MODERATE_LOW",sampOrPat='sample',prefix=''){
  #tab <- read.table(synfile@filePath,sep = '\t',header = T)
  
 # if(!includeSilent)
  #  tab <- tab[which(tab$Mutation_Type!='Silent'),]
  tab$Patient=sapply(tab$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
  #mut.idx = which(sapply(tab$Mutation_Status,tolower)%in%mutationType)
  tab <- subset(tab,Effect%in%unlist(strsplit(effect,split='_')))
  tab<- subset(tab,Status%in%mutationType)
  
 # tab<- tab[which(tab$Hugo_Symbol%in%expressed.genes),]
  if(sampOrPat=='sample')
    num.patients = tab %>% group_by(Gene) %>% summarize(Samples = n_distinct(Sample))
  else
    num.patients = tab%>% group_by(Gene)%>%summarize(Patients=n_distinct(Patient))
  
  fname= paste(prefix,paste(mutationType,collapse='_and_'),'mutationsWith',effect,'effectNum',sampOrPat,'perGene.tab',sep='')  
  write.table(num.patients,file=fname,col.names=F,row.names=F,quote=F)
  return(fname)
}

files=c()
for(mutType in list(c('Germline'),c("StrongSomatic","LikelySomatic"),c("StrongSomatic","LikelySomatic","StrongLOH","LikelyLOH"))){
#mutType<-'germline'
 for(mt in c("HIGH","HIGH_MODERATE","HIGH_MODERATE_LOW")){
      files=c(files,getCountsFromTable(expr.gene.muts05,mutationType=mutType,effect=mt,sampOrPat=ifelse(length(mutType)==1,'patient','sample'),prefix='p05'))
      files=c(files,getCountsFromTable(expr.gene.muts1,mutationType=mutType,effect=mt,sampOrPat=ifelse(length(mutType)==1,'patient','sample'),prefix='p1'))
      
       }
  
}

this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/bin/formatHotnetDataFromVardict.R'
for(f in files){
  synStore(File(f,parentId='syn6128017'),executed=list(list(url=this.script)))
}
