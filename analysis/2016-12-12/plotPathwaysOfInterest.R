###plot scores
source("../../bin/wgsAnalysis.R")
source("../../analysis/2016-09-20/geneSampleMatrix.R")
##takes a while, don't source if not plotting ssGSEAandGSVA
source("../../bin/ssGSEAandGSVA.R")
library(ggplot2)

p05='syn6097853'
p1='syn6099307'

this.file='https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2016-12-12/plotPathwaysOfInterest.R'

plotPathwayAcrossScores<-function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes[-1]
  if(TRUE %in% index) {
  pathwayMutated <- apply(patient.sample.muts[index, overlap, drop = FALSE], 2, any)
  ##add patient to data frame
  df<-data.frame(pathwayMutated=pathwayMutated ,score=patient.sample.vars[overlap],patient=sapply(overlap,function(x) paste(unlist(strsplit(x,split=' '))[1:2],collapse=' ')))
  pval<-'1.0'
  ##calculate p-value
  try(pval<-format(t.test(score~pathwayMutated,data=df)$p.value,digits=4), silent = TRUE)
  ##make color by patient and shape by mutation
  p<-ggplot()+geom_boxplot(data=df,aes(x=pathwayMutated,y=score), outlier.color = "NA", color = pathwayMutated)+
    scale_color_manual(c(TRUE=="green", FALSE=="blue"))+
    geom_jitter(data=df,aes(x=pathwayMutated,y=score,shape=pathwayMutated,color=patient))+
    ggtitle(paste(pathwayGenes[1],collapse="_",'mutation status by\n',testName,'scores\np=',pval))+
    theme(plot.title=element_text(hjust=0.5))
  print(p)
  #p<-ggplot()+geom_boxplot(data=df,aes(x=pathwayMutated,y=score))+geom_jitter(data=df,aes(x=pathwayMutated,y=score,color=pathwayMutated))+ggtitle(paste(pathwayGenes[1],collapse="_"),'mutation status by\n',testName,'scores')
  ggsave(paste(paste(pathwayGenes[1],collapse="_"),'mutationsBy',testName,'score.png',sep='_'))
  file.name <- paste(pathwayGenes[1],'mutationsBy',testName,'score.png',sep='_')
  #synStore(File(file.name, parentId='syn7874685'), executed=this.file)
  return(df)
  } else {
    print(paste(pathwayGenes[1], " not in seq data"))
  }
}


calculatePvalues <- function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName) {
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to evaluate significant difference'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes
  pathwayMutated <- apply(patient.sample.muts[index, overlap, drop = FALSE], 2, any)
  df<-data.frame(pathwayMutated = pathwayMutated, score = patient.sample.vars[overlap])
  pval<-1
  try((pval <- t.test(df$score ~ df$pathwayMutated)$p.value), silent = TRUE)
  print(pval)
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

pathway.list <- data.frame("CREBBP" = "CREBBP")

pvalGermline.df <- data.frame("Signature"=c("placeholder"), "Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalSomatic.df <- data.frame("Signature"=c("placeholder"), "Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalEstimateGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalEstimateSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

##pdf('plotPathwaysAcrossMutScores.pdf')
for(x in pathway.list){
  for(cp in cib.paths){
    cs<-cib.scores[,cp]
    names(cs)<-rownames(cib.scores)
    #plotPathwayAcrossScores(cs,gl.vars.by.sample,x,paste(cp,'Germline',sep=''))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    signaturename <- print(cp[1], max.levels = 0)
    signaturename <- sub(pattern = "[1]", signaturename, replacement = "", fixed = TRUE)
    pvalGermline <- c(signaturename, pathwayname, calculatePvalues(cs,gl.vars.by.sample,x,paste('Germline',sep='')))
    names(pvalGermline) <- c("Signature", "Pathway", "p.value")
    pvalGermline.df <- rbind(pvalGermline.df, pvalGermline)
  ##########
    #plotPathwayAcrossScores(cs,som.vars,x,paste(cp,'Somatic',sep='')) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    signaturename <- print(cp[1], max.levels = 0)
    signaturename <- sub(pattern = "[1]", signaturename, replacement = "", fixed = TRUE)
    pvalSomatic <- c(signaturename, pathwayname, calculatePvalues(cs,som.vars,x,paste('Somatic',sep='')))
    names(pvalSomatic) <- c("Signature", "Pathway", "p.value")
    pvalSomatic.df <- rbind(pvalSomatic.df, pvalSomatic)
  } }

for(x in pathway.list){
  es<-as.numeric(est.scores['ESTIMATEScore',])
  names(es)<-colnames(est.scores)
  ##########
    plotPathwayAcrossScores(es,gl.vars.by.sample,x,paste('EstimateGermline',sep=''))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalEstimateGermline <- c(pathwayname, calculatePvalues(es,gl.vars.by.sample,x,paste('EstimateGermline',sep='')))
    names(pvalEstimateGermline) <- c("Pathway", "p.value")
    pvalEstimateGermline.df <- rbind(pvalEstimateGermline.df, pvalEstimateGermline)
  #########
    plotPathwayAcrossScores(es,som.vars,x,paste('EstimateSomatic',sep='')) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalEstimateSomatic <- c(pathwayname, calculatePvalues(es,som.vars,x,paste('EstimateSomatic',sep='')))
    names(pvalEstimateSomatic) <- c("Pathway", "p.value")
    pvalEstimateSomatic.df <- rbind(pvalEstimateSomatic.df, pvalEstimateSomatic)
}
   ##dev.off()

for(x in pathway.list){
  for(y in rownames(hallmark.ssGSEA)){
  hall<-as.numeric(hallmark.ssGSEA[y,])
  names(hall)<-colnames(hallmark.ssGSEA)
  ##########
  plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y))
  pathwayname <- capture.output(print(x[1], max.levels = 0))
  pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
  pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y)))
  names(pvalhallmarkGermline) <- c("Pathway", "p.value")
  pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
  #########
  plotPathwayAcrossScores(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)) 
  pathwayname <- capture.output(print(x[1], max.levels = 0))
  pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
  pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)))
  names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
  pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
}
}
pvalGermline.df <- pvalGermline.df[-1,]
pvalSomatic.df <- pvalSomatic.df[-1,]
pvalEstimateGermline.df<- pvalEstimateGermline.df[-1,]
pvalEstimateSomatic.df <- pvalEstimateSomatic.df[-1,]
pvalhallmarkGermline.df <- pvalhallmarkGermline.df[-1,]
pvalhallmarkSomatic.df <- pvalhallmarkSomatic.df[-1,]

pvalGermline.BH<- dplyr::filter(pvalGermline.df, p.value<1)
pvalGermline.BH$BH <- p.adjust(pvalGermline.BH$p.value, method = "BH")

pvalSomatic.BH<- dplyr::filter(pvalSomatic.df, p.value<1)
pvalSomatic.BH$BH <- p.adjust(pvalSomatic.BH$p.value, method = "BH")

pvalEstimateGermline.BH<- dplyr::filter(pvalEstimateGermline.df, p.value<1)
pvalEstimateGermline.BH$BH <- p.adjust(pvalEstimateGermline.BH$p.value, method = "BH")

pvalEstimateSomatic.BH<- dplyr::filter(pvalEstimateSomatic.df, p.value<1)
pvalEstimateSomatic.BH$BH <- p.adjust(pvalEstimateSomatic.BH$p.value, method = "BH")

pvalhallmarkGermline.BH <- dplyr::filter(pvalhallmarkGermline.df, p.value<1)
pvalhallmarkGermline.BH$BH <- p.adjust(pvalhallmarkGermline.BH$p.value, method = "BH")

pvalhallmarkSomatic.BH <- dplyr::filter(pvalhallmarkSomatic.df, p.value<1)
pvalhallmarkSomatic.BH$BH <- p.adjust(pvalhallmarkSomatic.BH$p.value, method = "BH")


write.table(pvalGermline.BH, file="pvalGermline.txt")
write.table(pvalSomatic.BH, file="pvalSomatic.txt")
write.table(pvalEstimateGermline.BH, file="pvalEstimateGermline.txt")
write.table(pvalEstimateSomatic.BH, file="pvalEstimateSomatic.txt")

pvalhallmarkSomatic.BH <- pvalhallmarkSomatic.BH[,order('BH')]

#synStore(File('pvalGermline.txt', parentId='syn7874685'), executed=this.file)
#synStore(File('pvalSomatic.txt', parentId='syn7874685'), executed=this.file)
#synStore(File('pvalEstimateGermline.txt', parentId='syn7874685'), executed=this.file)
#synStore(File('pvalEstimateSomatic.txt', parentId='syn7874685'), executed=this.file)

