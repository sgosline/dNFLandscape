#recalculate mutation files to include all data
source("../../../dermalNF/bin/WGSData_VarDict.R")
library(parallel)

res<-mclapply(c('1','2','3','4','5','6','8','9','11'),function(patientId){
    r1=storeMutsForAllGenes(impact=c("HIGH","LOW","MODERATE","","MODIFIER"),pval=0.1,patientNumber=patientId)
    r2=storeMutsForAllGenes(impact=c("HIGH","LOW","MODERATE","","MODIFIER"),pval=0.05,patietnNumber=patientId)
},mc.cores=6)
