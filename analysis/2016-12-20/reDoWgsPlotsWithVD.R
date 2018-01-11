source("../../bin/wgsAnalysis.R")
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-20/reDoWgsPlotsWithVD.R'
p05='syn6097853'
p1='syn6099307'
##now try a workflow


#visualize in heatmap
allfiles<-c()

for(vd in c(10,20,50)){
high01<-subset(expr.gene.muts1,Effect=='HIGH'&VD>vd)
high05<-subset(expr.gene.muts05,Effect=='HIGH'&VD>vd)
res=plotMutsAcrossSamples(subset(high01,Status=='Germline'&Gene%in%cancer.genes$Gene.Symbol),samples=F,minVal=1,prefix=paste('GLmutsCosmicHighEffect_VDgt',vd,sep=''))
res.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=0,prefix=paste('SommutsCosmicP1HighEffect_VDgt',vd,sep=''))
res.s.05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")&Gene%in%cancer.genes$Gene.Symbol),minVal=1,prefix=paste('SommutsCosmicP05HighEffect_VDgt',vd,sep=''))
allfiles<-c(res$file)#,res.s$file,res.s.05$file)

res.all=plotMutsAcrossSamples(subset(high01,Status=='Germline'),samples=F,minVal=6,prefix=paste('GLmutsHighEffect_VDgt',vd,sep=''))
res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix=paste('SommutsP1HighEffect_VDgt',vd,sep=''))
res.all.s05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=2,prefix=paste('SommutsP05HighEffect_VDgt',vd,sep=''))
allfiles<-c(allfiles,res.all$file,res.all.s$file,res.all.s05$file)

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=4,prefix=paste('SommutsP1HighEffect_VDgt',vd,sep=''))
res.all.s05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=4,prefix=paste('SommutsP05HighEffect_VDgt',vd,sep=''))
allfiles<-c(allfiles,res.all.s$file,res.all.s05$file)

res.all.s=plotMutsAcrossSamples(subset(high01,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix=paste('SommutsP1HighEffect_VDgt',vd,sep=''))
res.all.s05=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix=paste('SommutsP05HighEffect_VDgt',vd,sep=''))
allfiles<-c(allfiles,res.all.s$file,res.all.s05$file)
}
#res05.all.s=plotMutsAcrossSamples(subset(high05,Status%in%c("StrongSomatic","LikelySomatic")),minVal=5,prefix='SommutsP05HighEffect')
for(file in allfiles){
  if(length(grep('p05',file)>0))
    mf<-p05
  else
    mf<-p1
  synStore(File(file,parentId='syn5605256'),used=list(list(entity=mf)),executed=list(list(url=this.script)))
}
