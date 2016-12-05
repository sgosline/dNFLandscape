#encode data processing

require(synapseClient)
synapseLogin()

gene.mapping<-read.table('../data/ensemblHugoGeneMapping.txt',header=T,sep='\t')
metadata<-read.table(synGet("syn6023670")@filePath,sep='\t',header=T)

##match the sample names to the metadata
getSampleNamesForMatrix<-function(matrix){
  inds=match(colnames(matrix),metadata$File.accession)
  df<-data.frame(Sample=metadata$Biosample.term.name[inds],
                 Type=metadata$Biosample.type[inds],
                 Library=metadata$Library.made.from[inds],
                 Age=metadata$Biosample.Age[inds],
                 Sex=metadata$Biosample.sex[inds],
                 Fraction=metadata$Biosample.subcellular.fraction.term.name[inds],
                 Replicate=metadata$Biological.replicate.s[inds])
  rownames(df)<-colnames(matrix)
  return(df)
                 
}

##match hugo gene names to the ensembl gene names
getGeneNamesForMatrix<-function(matrix){
  print('Getting gene names')
  inds=sapply(rownames(matrix),function(x) match(unlist(strsplit(x,split='.',fixed=T))[1],gene.mapping[,1]))
  gene.names<-as.character(gene.mapping[inds,2])
  navals<-which(is.na(gene.names))
  bvals=which(gene.names=='')
  
  rmatrix=matrix[-c(navals,bvals),]
  gn<-gene.names[-c(navals,bvals)]
  
  dups=which(duplicated(gn))
  rmatrix=rmatrix[-dups,]
  gn=gn[-dups]
  rownames(rmatrix)<-gn
  return(rmatrix)
  
}

#get the actual matrix
getEncodeSkinMatrix<-function(metric=c('TPM','FPKM'),alignment=c('hg19','grch38'),doVoomNorm=T){
  print(paste('Getting',metric,'of skin genes aligned to',alignment))
  if(metric=='TPM')
    if(alignment=='hg19')
      si='syn6035991'
    else
      si='syn6035993'
  else
    if(alignment=='hg19')
      si='syn6035999'
    else
      si='syn6035996'
  tab<-read.table(synGet(si)@filePath,sep=',',header=T)
  if(doVoomNorm){
          print("Performing VOOM normalization")
          library(limma)
          tab = voomWithQualityWeights(tab)$E
        }
  return(tab)

}

clusterSamples<-function(metric='TPM',alignment='grch38'){
  library(ggbiplot)
  matrix=getGeneNamesForMatrix(getEncodeSkinMatrix(metric,alignment))
  
  zv=which(apply(matrix,1,var)==0)
 # zc<-which(apply(matrix,2,var)==0)
  print('Clustering samples in matrix')
  if(length(zv)>0){
    matrix<-matrix[-zv,]
  }
  pc=prcomp(t(matrix),center=T,scale=T)
  samps<-getSampleNamesForMatrix(matrix)
  print(head(samps))
  png(paste('pcPlotsForEncodeSkin',metric,'valuesAlignedTo',alignment,'.png',sep=''))
  
  fname=paste('pcPlotsForEncodeSkin',metric,'valuesAlignedTo',alignment,'.png',sep='')
  ##plot by sample name
  p<-ggbiplot(pc,groups=samps$Sample,var.axes=F)
  ggsave(paste('bySample',fname,sep=''),p)
  
  #cell type
  p<-ggbiplot(pc,groups=samps$Type,var.axes=F)
  ggsave(paste('byCell',fname,sep=''),p)
  
  #library type
  p<-ggbiplot(pc,groups=samps$Library,var.axes=F)
  ggsave(paste('byLibrary',fname,sep=''),p)
  
  #fractions
  p<-ggbiplot(pc,groups=samps$Fraction,var.axes=F)
  ggsave(paste('byFraction',fname,sep=''),p)
  
  
}
