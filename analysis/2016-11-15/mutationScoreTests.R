source("../../bin/wgsAnalysis.R")
source("../../bin/geneSampleMatrix.R")
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-11-15/mutationScoreTests.R'
p05='syn6097853'
p1='syn6099307'

som.vars<-somaticGeneSampleMatrix()
germ.vars<-germlineGeneSampleMatrix()

#som.vars.df<-getMutsAcrossGenes(expr.gene.muts1,effect=c("HIGH"),germLine=c("Germline"),
#                                         som=c("LikelySomatic",'StrongSomatic'))

##re-plot somatic mutations with improved cutoffs
genes.with.more.vars<-names(which(apply(som.vars,1,function(x) length(which(x)))>2))
vars.per.gene<-rowSums(som.vars[genes.with.more.vars,]*1)
genes.with.more.vars<-genes.with.more.vars[order(vars.per.gene)]

samp.mut.burden<-colSums(som.vars*1)
genes.mutated<-rowSums(som.vars*1)
pheatmap(som.vars[genes.with.more.vars,]*1,cellwidth = 10,cellheight=10,cluster_cols = F,
          cluster_rows=F,annotation_col = data.frame(MutationalBurden=samp.mut.burden),
         annotation_row=data.frame(TimesMutated=genes.mutated),
         main='Samples and genes with >2 mutations',file='genesWithMoreThan2Muts.png')

##let's plot top 10?
full.stats<-sapply(rev(genes.with.more.vars)[1:10],function(gene){
  try(stats<-getMutationStatsForGene(expr.gene.muts1,gene=gene,doPlot=TRUE,effect=c("HIGH"),prefix='p1'))})


pngfiles=list.files('.')[grep('png',list.files('.'))]

for(pf in pngfiles){
  synStore(File(pf,parentId='syn7542129'),used=list(list(entity='syn6099307')),executed=list(list(url=this.script)))
}
