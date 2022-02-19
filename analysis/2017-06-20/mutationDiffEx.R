##do differential expression analysis on the samples that have specific germline/somatic variants of genes.


#source("../../bin/wgsAnalysis.R")
source("../2016-09-20/geneSampleMatrix.R")
source("../../../dermalNF/bin/dermalNFData.R")
library(limma)
##select based on individual mutations
#genes=c("CDC27","CREBBP")
genes=c("CREBBP")
##get gene expression
rna.counts<-rna_count_matrix(stored=TRUE,doNorm=FALSE,minCount=1,doLogNorm=FALSE,doVoomNorm=TRUE)
rna.fpkm<-rna_fpkm_matrix(byIsoform=FALSE)
colnames(rna.counts)<-tolower(patient_tumor_number_rna(colnames(rna.counts),quant='featureCounts'))
colnames(rna.fpkm)<-tolower(patient_tumor_number_rna(colnames(rna.fpkm),quant='cuffLinks'))

##get mutation data
som.vars<-somaticGeneSampleMatrix()

##select sample ids for which that gene is mutated
diffExByGene<-function(genelist,rna.mat){
  #create vector
  g=paste(genelist,collapse='_')

  print(paste("Separating samples by",g))
  mut.vec<-rep(paste('WT',g),ncol(rna.mat))
  names(mut.vec)<-colnames(rna.mat)
  
  if(length(genelist)==1)
    vars<-intersect(names(which(som.vars[g,])),colnames(rna.mat))
  else
    vars<-intersect(colnames(som.vars)[unique(unlist(apply(som.vars[genelist,],1,which)))],colnames(rna.mat))
  
  mut.vec[vars]<-rep(paste('MUT',g),length(vars))
  
  mut.fac<-factor(mut.vec)
  levs<-levels(mut.fac)
  design=model.matrix(~mut.fac)
  
  colnames(design)[1:2]=levs
  fit <- lmFit(rna.mat, design)
  fit <- eBayes(fit)
  tab <- topTable(fit, coef=levs[length(levs)],p.value=0.1)
  
  print(paste('Found',nrow(tab),'differentially expressed genes for',g))
  return(tab)
  
}

require(dplyr)

plotGenesByMut<-function(mut.genes,expr.genes,rna.mat){
  
  ##create mutant column
  #create vector
  g=paste(mut.genes,collapse='_')
  
  expr.genes<-intersect(expr.genes,rownames(rna.mat))
  print(paste("Separating samples by",g))
  mut.vec<-rep('WT',ncol(rna.mat))
  names(mut.vec)<-colnames(rna.mat)
  
  if(length(mut.genes)==1)
    vars<-intersect(names(which(som.vars[g,])),colnames(rna.mat))
  else
    vars<-intersect(colnames(som.vars)[unique(unlist(apply(som.vars[mut.genes,],1,which)))],colnames(rna.mat))
  
  mut.vec[vars]<-rep('MUT',length(vars))
  names(mut.vec)<-colnames(rna.mat)
  
  ##now collect gene expression values for same
#  gene.expr<-rna.mat[expr.genes,]
  
  aug.mat<-data.frame(Gene=rownames(rna.mat),rna.mat)
  df<-tidyr::gather(aug.mat,'Patient','RNA Counts',2:ncol(aug.mat))
  df$Patient<-sapply(df$Patient,function(x) gsub('.',' ', x,fixed=T))
  mdf<-data.frame(Patient=names(mut.vec),Status=mut.vec)
  
  fulldf<-left_join(df,mdf,by='Patient')
  fulldf$`RNA Counts`=fulldf$`RNA Counts`+1
  
  mindf<-subset(fulldf,Gene%in%expr.genes)
  ##now that we have the proper data frame we can plot
g<-ggplot(mindf)+geom_boxplot(aes(x=Gene,y=`RNA Counts`,fill=Status))+ggtitle(paste('Gene expression  of samples divided by mutations of\n',paste(mut.genes,collapse=' ')))
  print(g)
  
}

counts.diffex<-diffExByGene(genes,rna.counts)
fpkm.diffex<-diffExByGene(genes,rna.fpkm)

up.list<-c(genes,'NFATC4')
tg.counts.diffex<-diffExByGene(up.list,rna.counts)
tg.fpkm.diffex<-diffExByGene(up.list,rna.fpkm)


full.list<-c(genes,"SSB",'MED25','MED12','MED17','MED16','TRIP4','IRF3','HOXB7','SS18L1','GATAD1','ITGA5','FLT4','SHC3','NTRK3','FRS2','FGFR1','SH3GL1','PIAS3','TGS1','CUX1','NFATC4','ANAPC5')
full.counts.diffex<-diffExByGene(up.list,rna.counts)
full.fpkm.diffex<-diffExByGene(up.list,rna.fpkm)

complement.system=c("C3","C2",'C5','C6','C7','C9','C1QA',
                    'C1QB','C1QC','C1R','C1S','C4A','C4B',
                    'C8A','C8B','C8G','C1R','C1S','C3AR1',
                    'C5AR1','CD46','CD59','CFB','CFD','CFH',
                    'CFHR1','CFHR2','CFHR3','CFHR4','CFHR5','CFI','CFP',
                    'CR1','CR1L','CR2','ITGAM','ITGAX','ITGB2')
                    
require(ggplot2)

#plotGenesByMut(c("CDC27","CREBBP"),c("C3","HRCT1"),rna.counts)
#ggsave("CDC27_CREBBP_diffexgenes_muts.png")

#plotGenesByMut(c("CDC27","CREBBP"),complement.system,rna.counts)

#ggsave("CDC27_CREBBP_allcomps_muts.png")


plotGenesByMut(mut.genes=c("CREBBP"),complement.system,rna.counts)

ggsave("CREBBP_allcomps_muts.png")

