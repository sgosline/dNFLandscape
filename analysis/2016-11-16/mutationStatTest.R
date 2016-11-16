source("../../bin/wgsAnalysis.R")
source("../../bin/geneSampleMatrix.R")
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-11-16/mutationStatTests.R'
p05='syn6097853'
p1='syn6099307'
#compare germline/somatic mutations to various features

som.vars<-somaticGeneSampleMatrix()
germ.vars<-germlineGeneSampleMatrix()

##now can we check for germline correlates. this function needs to be
##expanded for additional statistical tests
checkForMutationCorrelations<-function(patient.sample.vars,patient.sample.muts,pvalthresh=0.1){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  tpvals<-apply(patient.sample.muts[,overlap],1,function(x) t.test(patient.sample.vars[overlap],x)$p.value)
  adj<-p.adjust(tpvals)
  
  return(names(which(adj<pvalthresh)))
  
}


##germline variants that correlate with somatic burden....
patients<-sapply(colnames(germ.vars),function(x) unlist(strsplit(x,split=' '))[2])

patient.vars<-sapply(unique(patients),function(x){
  apply(germ.vars[,which(patients==x)],1,function(y) any(y))
})

gl.vars.by.sample<-patient.vars[,patients]
colnames(gl.vars.by.sample)<-colnames(germ.vars)

##plot estimate scores by mutation status
est.scores<-read.table(synGet('syn5908274')@filePath)
##get patient scores
colnames(est.scores)<-tolower(patient_tumor_number_rna(colnames(est.scores),quant='featureCounts'))

##now iterate through every score, and every germline variant. 
germ.est.vals<-apply(est.scores,1,function(x) checkForMutationCorrelations(x,gl.vars.by.sample,0.001))
##now iterate through every score and every somatic variant
som.est.vals<-apply(est.scores,1,function(x) checkForMutationCorrelations(x,som.vars,0.001))

##plot cibersort scores by mutation status

cib.scores<-read.table(synGet('syn5809355')@filePath,sep=',',header=T)
rownames(cib.scores)<-tolower(patient_tumor_number_rna(sapply(cib.scores$Input.Sample,function(x) gsub('”','',x))))



##do more correlation tests