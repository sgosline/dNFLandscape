##replot cibersort
require(synapseClient)
synapseLogin()

cs.res<-read.table(synGet('syn5809355')@filePath,sep=',',header=T)

##the PBK ids are missing from table, so need to query annotations
res<-synQuery("select patientID,tissueID,sampleID from entity where parentId=='syn5492805'")
#map<-unique(res)

#from table get generic tumor id
tres<-synTableQuery("SELECT Patient,RnaID,TumorNumber,'RNASeq (Cufflinks)' FROM syn5556216 where RnaID is not NULL")@values

idx<-match(res$entity.id,tres$`RNASeq (Cufflinks)`)

dres<-res[which(!is.na(idx)),]
tres<-tres[idx[which(!is.na(idx))],]

full.map<-cbind(dres,tres)

#map tumors to sample ids
sampleIds<-sapply(cs.res$Input.Sample,function(x){
  y=which(full.map$entity.sampleID==gsub('”','',x))
  paste("Patient",full.map$Patient[y],"Tumor",full.map$TumorNumber[y])
  
})

rownames(cs.res)<-sampleIds
cs.res<-cs.res[,-1]

this.file='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-09-16/rePlotCiberSort.R'
library(pheatmap)
pheatmap(cs.res[,1:22],annotation_row=cs.res[,23:25],cellheight=10,cellwidth = 10,filename='ciberSortRePlotted.png')
synStore(File('ciberSortRePlotted.png',parentId='syn5809348'),executed=list(list(url=this.file)))


