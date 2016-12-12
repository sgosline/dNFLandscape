
###plot scores
source("../../bin/wgsAnalysis.R")
source("../../bin/geneSampleMatrix.R")
library(ggplot2)
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-02/plotGenesOfInterest.R'
p05='syn6097853'
p1='syn6099307'
#compare germline/somatic mutations to various f

plotGeneAcrossScores<-function(patient.sample.vars,patient.sample.muts,geneName,testName){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  if(!geneName%in%rownames(patient.sample.muts))
     return(data.frame())
  df<-data.frame(geneMutated=patient.sample.muts[geneName,overlap],score=patient.sample.vars[overlap])
  #colnames(df)<-c(paste(geneName,'mutation'),paste(testName,'score'))
  p<-ggplot()+geom_boxplot(data=df,aes(x=geneMutated,y=score))+geom_jitter(data=df,aes(x=geneMutated,y=score,color=geneMutated))+ggtitle(paste(geneName,'mutation status by\n',testName,'scores'))
  print(p)
  ggsave(paste(geneName,'mutationsBy',testName,'score.png',sep='_'))
  return(df)
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

####
genelist=c('RAD9A','APTX','APEX1','GEN1','HIBCH','TBCK','SOS1','RABGAP','SGSM3','DOCK4','MUC1','BRD7','STAG2','CTDP1','CDC27','CREBBP')
cib.paths=c("Mast.cells.resting","Macrophages.M2")
#est.paths=c("ESTIMATEScore")

for(gene in genelist){
  for(cp in cib.paths){
    cs<-cib.scores[,cp]
    names(cs)<-rownames(cib.scores)
    plotGeneAcrossScores(cs,gl.vars.by.sample,gene,paste(cp,'Germline',sep='')) 
    plotGeneAcrossScores(cs,som.vars,gene,paste(cp,'Somatic',sep='')) 
  }
  es<-as.numeric(est.scores['ESTIMATEScore',])
  names(es)<-colnames(est.scores)
  plotGeneAcrossScores(es,gl.vars.by.sample,gene,paste('EstimateGermline',sep='')) 
  plotGeneAcrossScores(es,som.vars,gene,paste('EstimateSomatic',sep='')) 
}
  


