#Test pathway enrichment of various sets of mutations
source("../../bin/wgsAnalysis.R")

#2#evaluate most prevalent mutations at high impact, somatic/germline
p05_mutcounts<-getMutsAcrossGenes(expr.gene.muts05,effect=c("HIGH","MODERATE"),germLine=c("Germline","Deletion"),
                                  som=c("Deletion","LikelyLOH","StrongLOH","LikelySomatic",'StrongSomatic'))

p1_mutcounts<-getMutsAcrossGenes(expr.gene.muts1,effect=c("HIGH","MODERATE"),germLine=c("Germline","Deletion"),
                                 som=c("Deletion","LikelyLOH","StrongLOH","LikelySomatic",'StrongSomatic'))

