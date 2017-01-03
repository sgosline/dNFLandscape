#get gene mutation matrix - all genes by all samples!!!
library(synapseClient)
synapseLogin()

library(dplyr)
if(!exists('fullfile'))
  fullfile<-read.table(synGet('syn6099307')@filePath,sep='\t',header=T)

#do one table query of tumor numbers
tum.nums<-synTableQuery(paste("SELECT TumorNumber,WGS FROM syn5556216 where WGS is not NULL"))@values
getTumorNumber<-function(sampname){
  allvals<-unlist(strsplit(sampname,split='_'))
  wgs=''
  tn<-''
  if(length(allvals)>2){
    wgs=allvals[4]
    tn<-tum.nums$TumorNumber[match(wgs,tum.nums$WGS)]
    if(is.na(tn))
      tn<-'NA'
  }
  return(gsub("_",' ',gsub(wgs,tn,sampname)))
}

somaticGeneSampleMatrix<-function(muts=fullfile,effect=c("HIGH"),vd=0){
    som.muts<-filter(muts,Status%in%c('LikelySomatic','StrongSomatic'))%>%filter(Effect%in%effect)
    samp.gene<-reshape2::acast(som.muts,Gene~Sample)
    colnames(samp.gene)<-sapply(colnames(samp.gene),getTumorNumber)

    bin.mat<-samp.gene>0
    return(bin.mat)
}

germlineGeneSampleMatrix<-function(muts=fullfile,effect=c("HIGH"),vd=0){
  gl.muts<-filter(muts,Status=='Germline')%>%filter(Effect%in%effect)
  samp.gene<-reshape2::acast(gl.muts,Gene~Sample)
  colnames(samp.gene)<-sapply(colnames(samp.gene),getTumorNumber)
  bin.mat<-samp.gene>0
  return(bin.mat)
}
