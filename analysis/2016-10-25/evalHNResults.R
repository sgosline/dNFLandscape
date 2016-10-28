source("../../bin/wgsAnalysis.R")

##can we plot individual genes to view
genes=c("MDC1","CDC27","NRP2","VEGFC","SEMA3B","FLT4","CREBBP")
genes=c("PRKDC",'MAP3K14','VPS11','RRM1','RECQL4',"CDC27","CREBBP")


for(g in genes){
  stats<-getMutationStatsForGene(expr.gene.muts1,gene=g,doPlot=TRUE,effect=c("HIGH"),prefix='p1')

}

##moderate impact mutations will be a tough sell.  can we find any high impact ones?  
#we should probably store some oft hose images on synapse. 
