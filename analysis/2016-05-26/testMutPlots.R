source("../../dermalNF/bin/WGSData_VarDict.R")
cancer.genes<-read.csv('../../data/Census_allTue Jan 19 18-58-56 2016.csv')
if(!exists('expr.gene.muts05'))
  expr.gene.muts05<-read.table(synGet("syn6097853")@filePath,sep='\t')

if(!exists('expr.gene.muts1'))
  expr.gene.muts1<-read.table(synGet("syn6099307")@filePath,sep='\t')

##for a particular gene, collect all the mutational details across all samples
getMutationStatsForGene<-function(expr.gene.muts,gene='NF1',doPlot=FALSE,effect=c("LOW","MODERATE","HIGH"),prefix='p05'){

  sdf<-subset(expr.gene.muts,Gene==gene)
 # if()
  sdf$Patient<-sapply(sdf$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
  sdf$Patient<-sapply(sdf$Patient,function(x) paste(x,'(n=',length(unique(sdf$Sample[which(sdf$Patient==x)])),')'))
  transcripts<-unique(sdf$Transcript)
  print(paste("Found",length(transcripts),'unique transcripts for',gene))
  sdf$Effect=factor(as.character(sdf$Effect),levels=c("LOW","MODERATE",'HIGH'))
  sdf<-subset(sdf,Effect%in%effect)

  
  if(nrow(sdf)==0)
    return(sdf)
  
  noaa=which(sdf$Amino_Acid_Change=='')
  if(length(noaa)>0){
    levels(sdf$Amino_Acid_Change)<-c(levels(sdf$Amino_Acid_Change),as.character(sdf$Codon_Change[noaa]))
    sdf$Amino_Acid_Change[noaa]<-sdf$Codon_Change[noaa]
  }
  fname=''
  if(doPlot){
    require(ggplot2)
    sdf$Amino_Acid_Change<-sapply(sdf$Amino_Acid_Change,function(x) if(nchar(as.character(x))>20) strtrim(as.character(x),20) else x)
    for(t1 in transcripts){
      tdf<-subset(sdf,Transcript==t1)
      #plot all
      p<-ggplot(tdf)+geom_jitter(aes(y=Patient,x=Amino_Acid_Change,col=Functional_Class,size=Effect))+facet_grid(.~Status)
 #   p  +facet_grid(.~Status)
      p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(paste('Gene',gene,'Transcript',t1))#+facet_grid(.~Status)
      fname=paste('gene',gene,'transcript',t1,'mutationsByType',prefix,'.png',sep='')
      ggsave(p,file=fname)
      
      #just germline/LOH/Deletion
      p<-ggplot(subset(tdf,Status%in%c('Germline')))+geom_point(aes(y=Patient,x=Amino_Acid_Change,col=Functional_Class,size=Effect))#+facet_grid(.~Status)
      #   p  +facet_grid(.~Status)
      p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(paste('Gene',gene,'Transcript',t1))#+facet_grid(.~Status)
      fname=paste('gene',gene,'transcript',t1,'GermlinemutationsByType',prefix,'.png',sep='')
      ggsave(p,file=fname)
      
      ##just somatic
      p<-ggplot(subset(tdf,Status%in%c("LikelySomatic",'StrongSomatic','LikelyLOH','StrongLOH','Deletion')))+geom_point(aes(y=Sample,x=Amino_Acid_Change,col=Functional_Class,size=Effect,shape=Status))#+facet_grid(.~Status)
      #   p  +facet_grid(.~Status)
      p<-p+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(paste('Gene',gene,'Transcript',t1))#+facet_grid(.~Status)
      fname=paste('gene',gene,'transcript',t1,'SomaticmutationsByType',prefix,'.png',sep='')
      ggsave(p,file=fname)
      
    }
  }
 return(list(data=sdf,file=fname)) 
}

##collect counts of all mutations, somatic and germline, across patients
getMutsAcrossGenes<-function(expr.gene.muts,effect=c("HIGH"),germLine=c("Germline"),
                               som=c("Deletion","LikelyLOH","StrongLOH","LikelySomatic",'StrongSomatic')){
    sdf<-subset(expr.gene.muts,Effect%in%effect)
    sdf$Patient<-sapply(sdf$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
    counts<-sdf%>%group_by(Gene,Status)%>%summarize(nPatients=n_distinct(Patient))%>%arrange(desc(nPatients))
    gl.counts<-subset(sdf,Status%in%germLine)%>%group_by(Gene)%>%summarize(nPatients=n_distinct(Patient))%>%arrange(desc(nPatients))
    so.counts<-subset(sdf,Status%in%som)%>%group_by(Gene)%>%summarize(nSamps=n_distinct(Sample),nPatients=n_distinct(Patient))%>%arrange(desc(nSamps))
  return(list(somatic=so.counts,germline=gl.counts,all=counts))    
}

plotMutsAcrossSamples<-function(sdf,samples=TRUE,minVal,prefix=''){
  sdf$Present=rep(1,nrow(sdf))
  sdf$Patient<-sapply(sdf$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
  if(samples){
    prefix=paste(prefix,'_samples',sep='')
    mat<-reshape2::acast(sdf,Sample~Gene,value.var="Present",fun.aggregate=mean,fill=0)
    }else{
    mat<-reshape2::acast(sdf,Patient~Gene,value.var="Present",fun.aggregate=mean,fill=0)
      prefix=paste(prefix,'_patients',sep='')
    }
  mat=mat[,which(colSums(mat)>=minVal)]
    fname=paste(prefix,'WithAtLeast',minVal,'mutations.png',sep='')
    pheatmap(mat,cluster_rows=T,cluster_cols=T,filename = fname, cellheight=10,cellwidth=10)
    return(list(genes=rownames(mat),file=fname))
}


##goal is to calculate mutational burden across samples, then return file names
mutationBurdenAcrossSamples<-function(expr.gene.muts,effect=c("HIGH"),prefix=''){
  require(ggplot2)
  sdf<-subset(expr.gene.muts,Effect%in%effect)
  sdf$Patient<-sapply(sdf$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
  sdf$newStatus<-factor(sapply(as.character(sdf$Status),function(x){
    if(x%in%c("Deletion","Germline")) {
      return("GermOrDel") 
    }else if(x%in%c("LikelyLOH","StrongLOH")){
      return("LOH")
      }else if(x%in%c("StrongSomatic","LikelySomatic")){
        return("Somatic")}
    else{return(x)}}))
  
  
  ##first classsify everything by original status
  genes_per_sample<-sdf%>% group_by(Sample,Status)%>%summarize(numGenes=n_distinct(Gene))%>%left_join(sdf,by='Sample')%>%select(Sample,numGenes,Patient,Status=Status.x)
  ggplot(unique(genes_per_sample))+geom_boxplot(aes(x=Patient,y=numGenes,fill=Status))
  f1=paste(prefix,'mutationalBurdenPerSampleByStatusWith',paste(effect,collapse='_'),'impact.png',sep='')
  ggsave(f1)
  ##then, acknowleding taht germline mutations can can be conflated with deletions, summarize
  genes_per_sample<-sdf%>% group_by(Sample,newStatus)%>%summarize(numGenes=n_distinct(Gene))%>%left_join(sdf,by='Sample')%>%select(Sample,numGenes,Patient,newStatus=newStatus.x)
  ggplot(unique(genes_per_sample))+geom_boxplot(aes(x=Patient,y=numGenes,fill=newStatus))
  f2=paste(prefix,'mutationalBurdenPerSampleByCombinedStatusWith',paste(effect,collapse='_'),'impact.png',sep='')
  ggsave(f2)  
  return(c(f1,f2))##return file names for upload to synapse
  
}

##goal is to calculate mutational burden across samples, then return file names
mutationBurdenAcrossPatients<-function(expr.gene.muts,effect=c("HIGH"),prefix=''){
  require(ggplot2)
  sdf<-subset(expr.gene.muts,Effect%in%effect)
  sdf$Patient<-sapply(sdf$Sample,function(x) paste(unlist(strsplit(as.character(x),split='_'))[1:2],collapse='_'))
  sdf$newStatus<-factor(sapply(as.character(sdf$Status),function(x) if(x%in%c("Deletion","Germline")) return("GermOrDel") else return(x)))
  
  ##first classsify everything by original status
  genes_per_patient<-sdf%>% group_by(Patient,Status)%>%summarize(numGenes=n_distinct(Gene))
  ggplot(unique(genes_per_patient))+geom_bar(aes(x=Patient,y=numGenes,fill=Status),stat='identity',position='dodge')
  f1=paste(prefix,'mutationalBurdenPerPatientByStatusWith',paste(effect,collapse='_'),'impact.png',sep='')
  ggsave(f1)
  ##then, acknowleding taht germline mutations can can be conflated with deletions, summarize
  genes_per_patient<-sdf%>% group_by(Patient,newStatus)%>%summarize(numGenes=n_distinct(Gene))
  ggplot(unique(genes_per_patient))+geom_bar(aes(x=Patient,y=numGenes,fill=newStatus),stat='identity',position='dodge')
  f2=paste(prefix,'mutationalBurdenPerPatientByCombinedStatusWith',paste(effect,collapse='_'),'impact.png',sep='')
  ggsave(f2)  
  
  return(c(f1,f2))##return file names for upload to synapse
  
}

##now try a workflow
#1#evaluate NF1 mutations at any impact
p05_nf1=getMutationStatsForGene(expr.gene.muts05,gene='NF1',doPlot=TRUE,effect=c("LOW","MODERATE","HIGH"),prefix='p05')
p1_nf1=getMutationStatsForGene(expr.gene.muts1,gene='NF1',doPlot=TRUE,effect=c("LOW","MODERATE","HIGH"),prefix='p1')

p05_files=c(p05_nf1$file)
p1_files=c(p1_nf1$file)

#2#evaluate most prevalent mutations at high impact, somatic/germline
p05_mutcounts<-getMutsAcrossGenes(expr.gene.muts05,effect=c("HIGH","MODERATE"),germLine=c("Germline","Deletion"),
                             som=c("Deletion","LikelyLOH","StrongLOH","LikelySomatic",'StrongSomatic'))
 
p1_mutcounts<-getMutsAcrossGenes(expr.gene.muts1,effect=c("HIGH","MODERATE"),germLine=c("Germline","Deletion"),
                                  som=c("Deletion","LikelyLOH","StrongLOH","LikelySomatic",'StrongSomatic'))

#visualize in heatmap
res=plotMutsAcrossSamples(subset(expr.gene.muts1,Status=='Germline'&Gene%in%cancer.genes$Gene.Symbol),samples=F,minVal=1,prefix='GLmutsCosmic')
res.s=plotMutsAcrossSamples(subset(expr.gene.muts1,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=1,prefix='SommutsCosmic')

res.all=plotMutsAcrossSamples(subset(expr.gene.muts1,Status=='Germline'),samples=F,minVal=6,prefix='GLmuts')
res.all.s=plotMutsAcrossSamples(subset(expr.gene.muts1,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='Sommuts')

p1_rabep=getMutationStatsForGene(expr.gene.muts1,gene='RABEP1',doPlot=TRUE,effect=c("HIGH","MODERATE","LOW"),prefix='p1')
#upload to synapse
p1_kmt2c=getMutationStatsForGene(expr.gene.muts1,gene='KMT2C',doPlot=TRUE,effect=c("HIGH","MODERATE","LOW"),prefix='p1')

#3#evaluate top somatic muts in dividually


##4##plot distribution of somatic/germline muts
p05SampFiles<-mutationBurdenAcrossSamples(expr.gene.muts05,prefix='p05')
p1SampFiles<-mutationBurdenAcrossSamples(expr.gene.muts1,prefix='p1')

p05PatFiles<-mutationBurdenAcrossPatients(expr.gene.muts05,prefix='p05')
p1PatFiles<-mutationBurdenAcrossPatients(expr.gene.muts1,prefix='p1')

p05_files=c(p05_files,p05SampFiles,p05PatFiles)
p1_files=c(p1_files,p1SampFiles,p1PatFiles)

#upload to synapse


