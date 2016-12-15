#recalculate mutation files to include all data
source("../../../dermalNF/bin/WGSData_VarDict.R")

r1=storeMutsForAllGenes(impact=c("HIGH","LOW","MODERATE","","MODIFIER"),pval=0.1)
r2=storeMutsForAllGenes(impact=c("HIGH","LOW","MODERATE","","MODIFIER"),pval=0.05)