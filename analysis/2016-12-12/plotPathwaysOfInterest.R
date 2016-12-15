###plot scores
source("../../bin/wgsAnalysis.R")
source("../../analysis/2016-09-20/geneSampleMatrix.R")
library(ggplot2)
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-02/plotGenesOfInterest.R'
p05='syn6097853'
p1='syn6099307'
#compare germline/somatic mutations to various f

pval.df <- data.frame()

plotPathwayAcrossScores<-function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes
  if(TRUE %in% index) {
  pathwayMutated <- apply(patient.sample.muts[index, overlap, drop = FALSE], 2, any)
  df<-data.frame(pathwayMutated = pathwayMutated, score = cs[overlap])
  pathwaylabel <- colnames(x)
  print(pathwaylabel)
  #colnames(df)<-c(paste(geneName,'mutation'),paste(testName,'score'))
  p<-ggplot()+geom_boxplot(data=df,aes(x=pathwayMutated,y=score))+geom_jitter(data=df,aes(x=pathwayMutated,y=score,color=pathwayMutated))+ggtitle(paste(paste(pathwayGenes,collapse="_"),'mutation status by\n',testName,'scores'))
  print(p)
  ggsave(paste(paste(pathwayGenes,collapse="_"),'mutationsBy',testName,'score.png',sep='_'))
  return(df)
  } else {
    print(paste(pathwayGenes, " not in seq data"))
  }
}

calculatePvalues <- function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName){

  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to evaluate significant difference'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes
  pathwayMutated <- apply(gl.vars.by.sample[index, overlap, drop = FALSE], 2, any)
  df<-data.frame(pathwayMutated = pathwayMutated, score = cs[overlap])
  pathwaylabel <- names(pathway.list)
  
#  if((sum(df[df==TRUE])>1) && (sum(df[df==FALSE])>1)) {
  pval <- t.test(df$score ~ df$pathwayMutated)$p.value
  pval <- c(pathwaylabel[x], pval)
  print(pval)
  pval.df <- rbind(pval.df, pval)
#  } else {
  print(paste(pathwayGenes, " lacking two mutation categories"))
}# }

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

pathway.list <- read.table('enrichr_mutated_GO_terms.txt', header = TRUE, sep = "\t")

pval.df <- data.frame()

pdf('plotPathwaysAcrossMutScores.pdf')
for(x in pathway.list){
  for(cp in cib.paths){
    cs<-cib.scores[,cp]
    names(cs)<-rownames(cib.scores)
    plotPathwayAcrossScores(cs,gl.vars.by.sample,x,paste(cp,'Germline',sep='')) 
    calculatePvalues(cs,gl.vars.by.sample,x,paste(cp,'Germline',sep=''))
    plotPathwayAcrossScores(cs,som.vars,x,paste(cp,'Somatic',sep='')) 
  }
  es<-as.numeric(est.scores['ESTIMATEScore',])
  names(es)<-colnames(est.scores)
  plotPathwayAcrossScores(es,gl.vars.by.sample,x,paste('EstimateGermline',sep='')) 
  plotPathwayAcrossScores(es,som.vars,x,paste('EstimateSomatic',sep='')) 
}
dev.off()


this.file='https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2016-12-12/plotPathwaysOfInterest.R'
synStore(File('plotPathwaysAcrossMutScores.pdf', parentId='syn7845166'), executed=this.file, used=list('syn5809355', 'syn5908274', p05, p1))







########scratch area####

print(paste('We have',length(overlap),'samples to check for mutation correlates'))
index <- row.names(gl.vars.by.sample) %in% pathway.list$GTPase_regulator_activity
if(TRUE %in% index) 
{
  pathwayMutated <- apply(gl.vars.by.sample[index, overlap, drop = FALSE], 2, any)
  df<-data.frame(pathwayMutated = pathwayMutated, score = cs[overlap])
  pathwaylabel <- names(pathway.list)
  print(pathwaylabel[1])
  
  #colnames(df)<-c(paste(geneName,'mutation'),paste(testName,'score'))
  p<-ggplot()+geom_boxplot(data=df,aes(x=pathwayMutated,y=score))+geom_jitter(data=df,aes(x=pathwayMutated,y=score,color=pathwayMutated))+ggtitle(paste(paste(pathwayGenes,collapse="_"),'mutation status by\n',testName,'scores'))
  print(p)
  ggsave(paste(paste(pathwayGenes,collapse="_"),'mutationsBy',testName,'score.png',sep='_'))
  if(any(df$pathwayMutated) && !all(df$pathwayMutated)) {
    pval <- t.test(df$score ~ df$pathwayMutated)$p.value
    pval.df <- rbind(pval.df, pathwayGenes = pval)
    return(df) }
  
  else {return(df)}
  
}
else{
  print(paste(pathwayGenes, " not in seq data"))
}
}

