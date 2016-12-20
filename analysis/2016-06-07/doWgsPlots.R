source("../../bin/wgsAnalysis.R")

##now try a workflow


#visualize in heatmap
high01<-subset(expr.gene.muts1,Effect=='HIGH')
high05<-subset(expr.gene.muts05,Effect=='HIGH')
res=plotMutsAcrossSamples(subset(high01,Status=='Germline'&Gene%in%cancer.genes$Gene.Symbol),samples=F,minVal=1,prefix='GLmutsCosmicHighEffect')
res.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=2,prefix='SommutsCosmicP1HighEffect')

res05.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=2,prefix='SommutsCosmicP05High05')

res.all=plotMutsAcrossSamples(subset(high01,Status=='Germline'),samples=F,minVal=6,prefix='GLmutsHighEffect')
res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP1HighEffect')
res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP05HighEffect')

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix='SommutsP1HighEffect')
res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix='SommutsP05HighEffect')

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=1,prefix='SommutsP1HighEffect')
res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=1,prefix='SommutsP05HighEffect')


high01<-subset(expr.gene.muts1,Effect%in%c("MODERATE",'HIGH'))
high05<-subset(expr.gene.muts05,Effect%in%c("MODERATE",'HIGH'))
res=plotMutsAcrossSamples(subset(high01,Status=='Germline'&Gene%in%cancer.genes$Gene.Symbol),samples=F,minVal=1,prefix='GLmutsCosmicHighModEffect')
res.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=2,prefix='SommutsCosmicP1HighModEffect')

res05.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=2,prefix='SommutsCosmicP05HighModEffect')

res.all=plotMutsAcrossSamples(subset(high01,Status=='Germline'),samples=F,minVal=6,prefix='GLmuts')
res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP1HighModEffect')
res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP05HighModEffect')

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix='SommutsP1HighModEffect')
res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix='SommutsP05HighModEffect')

for(g in c("MAML2","CREBBP","NOTCH2","RABEP1","SUZ12"))
  getMutationStatsForGene(expr.gene.muts1,gene=g,doPlot=T,prefix='p1')


