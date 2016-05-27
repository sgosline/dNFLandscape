source("../../bin/wgsAnalysis.R")

##now try a workflow

#1#evaluate NF1 mutations at any impact
p05_nf1=getMutationStatsForGene(expr.gene.muts05,gene='NF1',doPlot=TRUE,effect=c("LOW","MODERATE","HIGH"),prefix='p05')
p1_nf1=getMutationStatsForGene(expr.gene.muts1,gene='NF1',doPlot=TRUE,effect=c("LOW","MODERATE","HIGH"),prefix='p1')

##collect files to upload. 
p05_files=c(p05_nf1$file)
p1_files=c(p1_nf1$file)

#visualize in heatmap
res=plotMutsAcrossSamples(subset(expr.gene.muts1,Status=='Germline'&Gene%in%cancer.genes$Gene.Symbol),samples=F,minVal=1,prefix='GLmutsCosmic')
res.s=plotMutsAcrossSamples(subset(expr.gene.muts1,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=1,prefix='SommutsCosmicP1')

res05.s=plotMutsAcrossSamples(subset(expr.gene.muts05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=1,prefix='SommutsCosmicP05')

res.all=plotMutsAcrossSamples(subset(expr.gene.muts1,Status=='Germline'),samples=F,minVal=6,prefix='GLmuts')
res.all.s=plotMutsAcrossSamples(subset(expr.gene.muts1,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP1')
res05.all.s=plotMutsAcrossSamples(subset(expr.gene.muts05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP05')

p1_rabep=getMutationStatsForGene(expr.gene.muts1,gene='RABEP1',doPlot=TRUE,effect=c("HIGH","MODERATE","LOW"),prefix='p1')
#upload to synapse
p1_kmt2c=getMutationStatsForGene(expr.gene.muts1,gene='KMT2C',doPlot=TRUE,effect=c("HIGH","MODERATE","LOW"),prefix='p1')

#3#evaluate top somatic muts in dividually


##4##plot distribution of somatic/germline muts
p05SampFiles<-mutationBurdenAcrossSamples(expr.gene.muts05,prefix='p05')
p1SampFiles<-mutationBurdenAcrossSamples(expr.gene.muts1,prefix='p1')

p05PatFiles<-mutationBurdenAcrossPatients(expr.gene.muts05,prefix='p05')
p1PatFiles<-mutationBurdenAcrossPatients(expr.gene.muts1,prefix='p1')

p05_files=c(p05_files,p05SampFiles,p05PatFiles,res05.s$file,res05.all.s$file)
p1_files=c(p1_files,p1SampFiles,p1PatFiles,res$file,res.s$file,res.all$file,res.all.s$file)

##now upload the files to synapse

