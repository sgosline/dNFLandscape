##now focus on enrichment of mutations

source('../../dermalNF/bin/dermalNFData.R')
source("geneSampleMatrix.R")


library(mHG)

calcWilcoxonEnrichment<-function(mut.matrix,eigenGenes){
  rownames(eigenGenes)<-sapply(rownames(eigenGenes),function(x) tolower(gsub(" ",'_',x)))
  overlap<-intersect(rownames(eigenGenes),colnames(mut.matrix))
  mut.matrix<-mut.matrix[,overlap]
  ##we have such minimal overlap
  sp<-which(apply(mut.matrix,1,function(x) length(which(x)))>1)
  mut.matrix<-mut.matrix[sp,]
  
  eigenGenes<-eigenGenes[overlap,]
  gene.vals<-sapply(rownames(mut.matrix),function(x){
    evals<-sapply(colnames(eigenGenes),function(y){
      p<-mHG.test(mut.matrix[x,order(eigenGenes[,y])])$p.value
      return(p)
    })
  })
  return(gene.vals)
  }

this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-09-20/calcMutEnrichment.R'
mut.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-09-20/geneSampleMatrix.R'
library(ggbiplot)
library(pheatmap)
iterateOverClusters<-function(){
  ##first let's figure out which WGCNA run to use.
  quant=c("cuffLinks",'featureCounts')
  filter=c("Conn",'Sd')
  sign=c('signed','unsigned')
  
  res<-c()
  for(q in quant){
    for(f in filter){
      for(s in sign){
    ##now get the eigen genes
    efname<-paste(q,tolower(f),'filtered',tolower(s),'clusterEigenGenes.tab',sep='_')
    eeid=allfiles$entity.id[match(efname,allfiles$entity.name)]
    etab<-read.table(synGet(eeid)@filePath,header=T)
    
    rownames(etab)<-patient_tumor_number_rna(rownames(etab),q)
    title=paste('Eigengenes for',s,f,'filtered',q,'measurements')
    ggbiplot(prcomp(etab),groups=sapply(rownames(etab),function(x) paste(unlist(strsplit(x,split=' '))[1:2],collapse=' ')))+ggtitle(title)
    fn=paste(gsub(' ','_',title),'_PCA.png',sep='')
    ggsave(fn)
    synStore(File(fn,parentId='syn7256344'),executed=list(list(url=this.script)),used=list(list(entity=eeid)))
    
    sig.muts<-list(germline=c(),somatic=c())
    for(muts in c('somatic','germline')){
      ##now compute enrichment
      if(muts=='somatic')
        mut.matrix<-somaticGeneSampleMatrix()
      else
        mut.matrix<-germlineGeneSampleMatrix()
      
      pvals<-calcWilcoxonEnrichment(mut.matrix,etab)
      sig.muts[[muts]]<-apply(pvals,1,function(x) paste(names(which(x<0.05)),collapse=','))
      
      mtitle<-paste(muts,'Mutations in',title)
      fn=paste(gsub(' ','_',mtitle),'_EnrichmentPvals.tsv',sep='')
      write.table(pvals,file=fn,sep='\t')
      synStore(File(fn,parentId='syn7256344'),executed=list(list(url=this.script),list(url=mut.script)),used=list(list(entity=eeid),list(entity='syn6097853')))
      
      fn=paste(gsub(' ','_',mtitle),'_EnrichmentHeatmap.png',sep='')
      pheatmap(pvals,main=mtitle,cellwidth=10,cellheight=10,filename=fn)
      synStore(File(fn,parentId='syn7256344'),executed=list(list(url=this.script),list(url=mut.script)),used=list(list(entity=eeid),list(entity='syn6097853')))
      
      ##last create a table of all significantly enriched genes per cluster
      
    }
    df<-data.frame(sig.muts)
    fn=paste(gsub(' ','_',title),'_enrichedMuts.tsv',sep='')
    write.table(df,file=fn,sep='\t')
    synStore(File(fn,parentId='syn7256344'),executed=list(list(url=this.script),list(url=mut.script)),used=list(list(entity=eeid),list(entity='syn6097853')))
    
    
          }
    }
  }
}

iterateOverClusters()