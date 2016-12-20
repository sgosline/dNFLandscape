source("../../dermalNF/bin/WGSData_VarDict.R")
source("../../bin/wgsAnalysis.R")
cancer.genes<-read.csv('../../data/Census_allTue Jan 19 18-58-56 2016.csv')
if(!exists('expr.gene.muts05'))
  expr.gene.muts05<-read.table(synGet("syn6097853")@filePath,sep='\t')

if(!exists('expr.gene.muts1'))
  expr.gene.muts1<-read.table(synGet("syn6099307")@filePath,sep='\t')


##now try a workflow
#1#evaluate NF1 mutations at any impact
p05_nf1=getMutationStatsForGene(expr.gene.muts05,gene='NF1',doPlot=TRUE,effect=c("LOW","MODERATE","HIGH"),prefix='p05')
p1_nf1=getMutationStatsForGene(expr.gene.muts1,gene='NF1',doPlot=TRUE,effect=c("LOW","MODERATE","HIGH"),prefix='p1')

p05_files=c(p05_nf1$file)
p1_files=c(p1_nf1$file)


#upload to synapse


