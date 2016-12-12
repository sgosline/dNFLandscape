source("../../bin/wgsAnalysis.R")
source("../../bin/geneSampleMatrix.R")
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-11-30/mutationStatTests.R'
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
  tpvals<-apply(patient.sample.muts[,overlap],1,function(x){
    if(length(unique(x))==1)
      return(1.0)
    else
        return(t.test(patient.sample.vars[overlap],x)$p.value)})
  
  tpvals<-sort(tpvals,decreasing=F)
  adj<-p.adjust(tpvals)
  
  #adj<-sort(adj, decreasing=F)
  print(head(adj))
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
germ.est.vals<-sapply(rownames(est.scores),function(x) {
  es<-est.scores[x,]
  names(es)<-colnames(est.scores)
  checkForMutationCorrelations(es,gl.vars.by.sample,0.05)
  })
##now iterate through every score and every somatic variant
som.est.vals<-sapply(rownames(est.scores),function(x) {
  es<-est.scores[x,]
  names(es)<-colnames(est.scores)
  checkForMutationCorrelations(es,som.vars,0.05)
  })

sapply(names(germ.est.vals),function(x) 
  write(germ.est.vals[[x]],sep='\n',file=paste('germlineVariants',x,'EstimateScores0.0001.txt',sep='_')))

sapply(names(som.est.vals),function(x) 
  write(som.est.vals[[x]],sep='\n',file=paste('somaticVariants',x,'EstimateScores0.0001.txt',sep='_')))

##plot cibersort scores by mutation status

cib.scores<-read.table(synGet('syn5809355')@filePath,sep=',',header=T)
rownames(cib.scores)<-tolower(patient_tumor_number_rna(sapply(cib.scores$Input.Sample,function(x) gsub('”','',x))))

tr<-c(1,which(apply(cib.scores,2,var)<0.0001))
##now iterate through every score, and every germline variant. 
germ.cib.vals<-sapply(colnames(cib.scores)[-tr],function(x) {
  cs<-cib.scores[,x]
  names(cs)<-rownames(cib.scores)
  checkForMutationCorrelations(cs,gl.vars.by.sample,0.05)})

##now iterate through every score and every somatic variant
som.cib.vals<-sapply(colnames(cib.scores)[-tr],function(x){
  cs<-cib.scores[,x]
  names(cs)<-rownames(cib.scores)
 checkForMutationCorrelations(cs,som.vars,0.05)})

sapply(names(germ.cib.vals),function(x) 
  write(germ.cib.vals[[x]],sep='\n',file=paste('germlineVariants',x,'CibersortScores0.05.txt',sep='_')))

sapply(names(som.cib.vals),function(x) 
  write(som.cib.vals[[x]],sep='\n',file=paste('somaticVariants',x,'CibersortScores0.05.txt',sep='_')))

txtfiles=list.files('.')[grep('Cibersort',list.files('.'))]
for(file in txtfiles){
  synStore(File(file,parentId='syn7806889'),used=list(list(entity='syn5809355'),list(entity=p1)),executed=list(list(url=this.script)))
}

txtfiles=list.files('.')[grep('Estimate',list.files('.'))]
for(file in txtfiles){
  synStore(File(file,parentId='syn7806889'),used=list(list(entity='syn5908274'),list(entity=p1)),executed=list(list(url=this.script)))
}
#upload files to synapse