library(synapseClient)
synapseLogin()

##cluster RNA Seq with ratner data

#get pnf data
pnfData<-read.table(synGet('syn7124098')@filePath,sep='\t',header=T)
colnames(pnfData)<- gsub('..',' (',gsub('.clone.',' clone)',colnames(pnfData),fixed=T),fixed=T)
phenoData<-read.table(synGet('syn7139168')@filePath,sep='\t',header=T)

#get ratner data
ratnerData<-read.table(synGet('syn5950004')@filePath,sep='\t',header=T)
ratnerPheno<-read.table(synGet('syn5950620')@filePath,sep='\t',header=T)

comm.genes<-intersect(rownames(ratnerData),rownames(pnfData))
full.dat<-cbind(ratnerData[comm.genes,],pnfData[comm.genes,])

zscore<-function(x) {(x-mean(x,na.rm=T))/sd(x)}

z.dat<-apply(full.dat,2,zscore)
r.dat<-apply(full.dat,2,rank)
require(limma)
l.dat<-normalizeBetweenArrays(full.dat)
##cluster pnF segments with serra cells
##now get all the segments..

library(ggbiplot)



rat.cell<-ratnerPheno$cellType[match(colnames(ratnerData),ratnerPheno[,1])]
wal.cell<-paste('pNF',phenoData$Genotype[match(colnames(pnfData),rownames(phenoData))])

celltypes<-c(as.character(rat.cell),wal.cell)
names(celltypes)<-c(colnames(ratnerData),colnames(pnfData))

pc<-prcomp(t(r.dat))

ggbiplot(pc,var.axes=F,groups=celltypes)
