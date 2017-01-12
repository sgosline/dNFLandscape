source("../../bin/wgsAnalysis.R")
source("../../bin/geneSampleMatrix.R")
library(dplyr)
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
    else if(length(which(x))<2)
      return(1.0)
    else
      try(return(t.test(patient.sample.vars[overlap],x)$p.value))})
  
  tpvals<-sort(tpvals,decreasing=F)
  adj<-p.adjust(tpvals)
  
  #adj<-sort(adj, decreasing=F)
  print(head(adj))
  return(names(which(adj<pvalthresh)))
  
}

checkForMutationCorrelationsFisher<-function(patient.sample.vars,patient.sample.muts,pvalthresh=0.1){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  tpvals<-apply(patient.sample.muts[,overlap],1,function(x){
    if(length(unique(x))==1)
      return(1.0)
    else if(length(which(x))<2)
      return(1.0)
    else
      try(return(fisher.test(patient.sample.vars[overlap],x)$p.value))})
  
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


checkSampsAgainstSom<-function(samp.mat,samp.name,som.vars,pval=0.05){
  colnames(samp.mat)<-sapply(colnames(samp.mat),function(x) tolower(gsub('.',' ',x,fixed=T)))
  samp.samps<-sapply(rownames(samp.mat),function(x){
    svals<-samp.mat[x,]
    names(svals)<-colnames(samp.mat)
    checkForMutationCorrelations(svals,som.vars,pval)
    
  })
  files<-sapply(names(samp.samps),function(x){
    genes<-samp.samps[[x]]
    if(length(genes)>0){
      fname<-paste('somaticVariants',x,samp.name,'Scores',pval,'.txt',sep='_')
      write(genes,sep='\n',file=fname)
      return(fname)
    }
  })
  return(files)
}

##chi squared test for categorical variables
checkSampsAgainstGerm<-function(samp.mat,samp.name,germ.vars,pval=0.05){
  colnames(samp.mat)<-sapply(colnames(samp.mat),function(x) tolower(gsub('.',' ',x,fixed=T)))
  samp.samps<-sapply(rownames(samp.mat),function(x){
    svals<-samp.mat[x,]
    names(svals)<-colnames(samp.mat)
    checkForMutationCorrelationsFisher(svals,germ.vars,pval)
    
  })
  files<-sapply(names(samp.samps),function(x){
    genes<-samp.samps[[x]]
    if(length(genes)>0){
      fname<-paste('germlineVariants',x,samp.name,'Scores',pval,'.txt',sep='_')
      write(genes,sep='\n',file=fname)
      return(fname)
    }
  })
  return(files)
}


##get cutaneous NF sample information
sample.info<-synTableQuery('SELECT * FROM syn5556216')@values
patient.info<-synTableQuery('SELECT * FROM syn7342635')@values

##prep tumor variables for comparison to somatic variants
sample.info<-dplyr::full_join(sample.info,patient.info, by="Patient")
sample.info<-dplyr::filter(sample.info, !is.na(WGS))
sample.info<-dplyr::filter(sample.info, !is.na(TumorNumber))

sample.info.somatic<-dplyr::select(sample.info, Patient, TumorNumber, Length_in_mm, Gender, Age, Inherited, Pain, Itching, rapidly_growing_lesions)
sample.info.somatic$Gender[sample.info.somatic$Gender=='Female'] <- 1
sample.info.somatic$Gender[sample.info.somatic$Gender=='Male'] <- 0
sample.info.somatic[6:9] <- apply(sample.info.somatic[6:9], 2, function(x) {x[x!='FALSE'] <- 1; x[x=='FALSE'] <- 0; x})

sample.info.somatic<-t(dplyr::mutate(sample.info.somatic, paste("Patient.", Patient, ".Tumor.", TumorNumber, sep = "")))
samplenames<-unname(c(sample.info.somatic[10,]))
colnames(sample.info.somatic) <- samplenames
sample.info.somatic<-as.matrix(sample.info.somatic[3:9,])
class(sample.info.somatic) <- 'numeric'

checkSampsAgainstSom(sample.info.somatic, 'Tumor_Info', som.vars, pval=0.005)

##prep tumor variables for comparison to somatic variants
sample.info.germline<-dplyr::select(sample.info, Patient, TumorNumber, 
                                    Length_in_mm, Race, Gender, Age, Menopause, 
                                    Inherited, Pregnancies, Births, Pain, Itching, 
                                    NumberOfPlexiforms, Manifestations, genetic_test_confirmation_nf1, 
                                    cafe_au_lait_spots, armpit_freckles, groin_freckles, lisch_nodules, 
                                    neurofibroma, age_of_neurofibromas, plexiform_neurofibroma, location_plexiform,
                                    history_of_low_vitamin_d, during_pregnancy, increased_during_pregnancy, 
                                    grew_during_pregnancy, family_memeber_with_nf1, rapidly_growing_lesions)

sample.info.germline<-t(dplyr::mutate(sample.info.germline, paste("Patient.", Patient, ".Tumor.", TumorNumber, sep = "")))
samplenames<-unname(c(sample.info.germline[30,]))
colnames(sample.info.germline) <- samplenames
sample.info.germline<-as.matrix(sample.info.germline[3:29,])

checkSampsAgainstGerm(sample.info.germline, 'Tumor_Info', germ.vars, pval=0.005)


