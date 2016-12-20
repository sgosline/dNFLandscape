source("../../bin/wgsAnalysis.R")
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-19/reDoWgsPlots.R'
p05='syn6097853'
p1='syn6099307'
##now try a workflow


#visualize in heatmap
high01<-subset(expr.gene.muts1,Effect=='HIGH')
high05<-subset(expr.gene.muts05,Effect=='HIGH')
allfiles<-c()
res=plotMutsAcrossSamples(subset(high01,Status=='Germline'&Gene%in%cancer.genes$Gene.Symbol),samples=F,minVal=1,prefix='GLmutsCosmicHighEffect')
res.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=1,prefix='SommutsCosmicP1HighEffect')
res.s.05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=1,prefix='SommutsCosmicP05HighEffect')
allfiles<-c(res$file,res.s$file,res.s.05$file)

res.all=plotMutsAcrossSamples(subset(high01,Status=='Germline'),samples=F,minVal=6,prefix='GLmutsHighEffect')
res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix='SommutsP1HighEffect')
res.all.s05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix='SommutsP05HighEffect')
allfiles<-c(allfiles,res.all$file,res.all.s$file,res.all.s05$file)

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=4,prefix='SommutsP1HighEffect')
res.all.s05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=4,prefix='SommutsP05HighEffect')
allfiles<-c(allfiles,res.all.s$file,res.all.s05$file)

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP1HighEffect')
res.all.s05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP05HighEffect')
allfiles<-c(allfiles,res.all.s$file,res.all.s05$file)

#res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP05HighEffect')
for(file in allfiles){
  if(length(grep('p05',file)>0))
    mf<-p05
  else
    mf<-p1
  synStore(File(file,parentId='syn5605256'),used=list(list(entity=mf)),executed=list(list(url=this.script)))
}
