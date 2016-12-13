###plot scores
source("../../bin/wgsAnalysis.R")
source("../../analysis/2016-09-20/geneSampleMatrix.R")
library(ggplot2)
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-02/plotGenesOfInterest.R'
p05='syn6097853'
p1='syn6099307'
#compare germline/somatic mutations to various f


###THIS FUNCTION REQUIRES 2 OR MORE GENES!!
plotPathwayAcrossScores<-function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes
  if(TRUE %in% index) 
  {
  pathwayMutated <- apply(patient.sample.muts[index, overlap, drop = FALSE], 2, any)
  df<-data.frame(pathwayMutated = pathwayMutated, score = cs[overlap])

  #colnames(df)<-c(paste(geneName,'mutation'),paste(testName,'score'))
  p<-ggplot()+geom_boxplot(data=df,aes(x=pathwayMutated,y=score))+geom_jitter(data=df,aes(x=pathwayMutated,y=score,color=pathwayMutated))+ggtitle(paste(paste(pathwayGenes,collapse="_"),'mutation status by\n',testName,'scores'))
  print(p)
  ggsave(paste(paste(pathwayGenes,collapse="_"),'mutationsBy',testName,'score.png',sep='_'))
  
  return(df) 
  }
  else{
    print(paste(pathwayGenes, " not in seq data"))
  }
}

som.vars<-somaticGeneSampleMatrix()
germ.vars<-germlineGeneSampleMatrix()

##germline variants that correlate with somatic burden....
patients<-sapply(colnames(germ.vars),function(x) unlist(strsplit(x,split=' '))[2])

patient.vars<-sapply(unique(patients),function(x){
  apply(germ.vars[,which(patients==x)],1,function(y) any(y))
})

gl.vars.by.sample<-patient.vars[,patients]
colnames(gl.vars.by.sample)<-colnames(germ.vars)

cib.scores<-read.table(synGet('syn5809355')@filePath,sep=',',header=T)
rownames(cib.scores)<-tolower(patient_tumor_number_rna(sapply(cib.scores$Input.Sample,function(x) gsub('”','',x))))

est.scores<-read.table(synGet('syn5908274')@filePath)
##get patient scores
colnames(est.scores)<-tolower(patient_tumor_number_rna(colnames(est.scores),quant='featureCounts'))

cib.paths=c("Mast.cells.resting","Macrophages.M2")

pathway.list <- na.omit(read.table('enrichr_mutated_GO_terms.txt', header = TRUE, sep = "\t"))

pdf('pathways.pdf')
for(x in pathway.list){
  for(cp in cib.paths){
    cs<-cib.scores[,cp]
    names(cs)<-rownames(cib.scores)
    plotPathwayAcrossScores(cs,gl.vars.by.sample,x,paste(cp,'Germline',sep='')) 
    plotPathwayAcrossScores(cs,som.vars,x,paste(cp,'Somatic',sep='')) 
  }
  es<-as.numeric(est.scores['ESTIMATEScore',])
  names(es)<-colnames(est.scores)
  plotPathwayAcrossScores(es,gl.vars.by.sample,x,paste('EstimateGermline',sep='')) 
  plotPathwayAcrossScores(es,som.vars,x,paste('EstimateSomatic',sep='')) 
}
dev.off()
