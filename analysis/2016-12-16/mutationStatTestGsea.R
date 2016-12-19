source("../../bin/wgsAnalysis.R")
source("../../bin/geneSampleMatrix.R")
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-16/mutationStatTestGsea.R'
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
    else if(which(x)<2)
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

gsva.hallmarks<-read.table(synGet('syn7862153')@filePath)
gsea.hallmarks<-read.table(synGet('syn7862155')@filePath)

gsva.immune<-read.table(synGet('syn7863864')@filePath)
gsea.immune<-read.table(synGet('syn7863609')@filePath)

gsva.onco<-read.table(synGet('syn7862147')@filePath)
gsea.onco<-read.table(synGet('syn7862149')@filePath)


checkSigsAgainstGerm<-function(sig.mat,sig.name,gl.vars.by.sample,pval=0.05){
  colnames(sig.mat)<-sapply(colnames(sig.mat),function(x) tolower(gsub('.',' ',x,fixed=T)))
  sig.sigs<-sapply(rownames(sig.mat),function(x){
    svals<-sig.mat[x,]
    names(svals)<-colnames(sig.mat)
    checkForMutationCorrelations(svals,gl.vars.by.sample,pval)
    
  })
  files<-sapply(names(sig.sigs),function(x){
    genes<-sig.sigs[[x]]
    if(length(genes)>0){
        fname<-paste('germlineVariants',x,sig.name,'Scores',pval,'.txt',sep='_')
        write(genes,sep='\n',file=fname)
        return(fname)
    }
  })
  return(files)
}

checkSigsAgainstSom<-function(sig.mat,sig.name,som.vars,pval=0.05){
  colnames(sig.mat)<-sapply(colnames(sig.mat),function(x) tolower(gsub('.',' ',x,fixed=T)))
  sig.sigs<-sapply(rownames(sig.mat),function(x){
    svals<-sig.mat[x,]
    names(svals)<-colnames(sig.mat)
    checkForMutationCorrelations(svals,som.vars,pval)
    
  })
  files<-sapply(names(sig.sigs),function(x){
    genes<-sig.sigs[[x]]
    if(length(genes)>0){
      fname<-paste('somaticVariants',x,sig.name,'Scores',pval,'.txt',sep='_')
      write(genes,sep='\n',file=fname)
      return(fname)
    }
  })
  return(files)
}

res<-checkSigsAgainstSom(gsva.hallmarks,'GSVA_Hallmarks',som.vars,pval=0.05)
res<-checkSigsAgainstSom(gsva.immune,'GSVA_Immune',som.vars,pval=0.05)
res<-checkSigsAgainstSom(gsva.onco,'GSVA_Onco',som.vars,pval=0.05)

res<-checkSigsAgainstSom(gsea.hallmarks,'GSEA_Hallmarks',som.vars,pval=0.05)
res<-checkSigsAgainstSom(gsea.immune,'GSEA_Immune',som.vars,pval=0.05)
res<-checkSigsAgainstSom(gsea.onco,'GSEA_Onco',som.vars,pval=0.05)
#upload files to synapse

res<-checkSigsAgainstSom(gsva.hallmarks,'GSVA_Hallmarks',som.vars,pval=0.01)
res<-checkSigsAgainstSom(gsva.immune,'GSVA_Immune',som.vars,pval=0.01)
res<-checkSigsAgainstSom(gsva.onco,'GSVA_Onco',som.vars,pval=0.01)

res<-checkSigsAgainstSom(gsea.hallmarks,'GSEA_Hallmarks',som.vars,pval=0.01)
res<-checkSigsAgainstSom(gsea.immune,'GSEA_Immune',som.vars,pval=0.01)
res<-checkSigsAgainstSom(gsea.onco,'GSEA_Onco',som.vars,pval=0.01)
#upload files to synapse