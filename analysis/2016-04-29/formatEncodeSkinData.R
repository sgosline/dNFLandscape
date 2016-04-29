##create annotated matrix from estimated skin counts
library(synapseClient)
synapseLogin()

allfiles<-synQuery("select name,id from entity where parentId=='syn6015781'")
metadata=read.table(synGet(allfiles$entity.id[grep('metadata',allfiles$entity.name)])@filePath,sep='\t',header=T,as.is=T)

gene.quants<-metadata[intersect(grep('tsv',metadata$File.format),grep('gene',metadata$Output.type)),]

##just get all the files, then annotate them and build into linear model later
tsvs<-gene.quants$File.accession
synids<-allfiles$entity.id[match(paste(tsvs,'tsv',sep='.'),allfiles$entity.name)]

all.tabs<-lapply(synids,function(x) {
  print(x)
  f<-read.table(synGet(x)@filePath,header=T)
  if(ncol(f)>4)
    return(f[,c('gene_id','TPM','FPKM')])
  else
    return(NULL)})

nv<-which(unlist(lapply(all.tabs,is.null)))
all.tabs<-all.tabs[-nv]
tsvs<-tsvs[-nv]
hg19.tsvs<-which(metadata[match(tsvs,metadata$File.accession),ncol(metadata)-1]=='hg19')
grch38.tsvs<-which(metadata[match(tsvs,metadata$File.accession),ncol(metadata)-1]=='GRCh38')

hg19.tpm.mat<-sapply(all.tabs[hg19.tsvs],function(x){
  counts<-as.numeric(x$TPM)
  names(counts)<-as.character(x$gene_id)
  return(counts)
})

hg19.fpkm.mat<-sapply(all.tabs[hg19.tsvs],function(x){
  counts<-as.numeric(x$FPKM)
  names(counts)<-as.character(x$gene_id)
  return(counts)
})
colnames(hg19.tpm.mat)<-tsvs[hg19.tsvs]
colnames(hg19.fpkm.mat)<-tsvs[hg19.tsvs]



grch38.tpm.mat<-sapply(all.tabs[grch38.tsvs],function(x){
  counts<-as.numeric(x$TPM)
  names(counts)<-as.character(x$gene_id)
  return(counts)
})

grch38.fpkm.mat<-sapply(all.tabs[grch38.tsvs],function(x){
  counts<-as.numeric(x$FPKM)
  names(counts)<-as.character(x$gene_id)
  return(counts)
})
colnames(grch38.tpm.mat)<-tsvs[grch38.tsvs]
colnames(grch38.fpkm.mat)<-tsvs[grch38.tsvs]

write.table(hg19.tpm.mat,file='hg19AlignedTPMs.csv',sep=',',quote=F)
write.table(hg19.fpkm.mat,file='hg19AlignedFPKMs.csv',sep=',',quote=F)

write.table(grch38.tpm.mat,file='grch38AlignedTPMs.csv',sep=',',quote=F)
write.table(grch38.fpkm.mat,file='grch38AlignedFPKMs.csv',sep=',',quote=F)

synStore(File('hg19AlignedTPMs.csv',parentId=''),used=list(url='',wasExecuted=TRUE))


