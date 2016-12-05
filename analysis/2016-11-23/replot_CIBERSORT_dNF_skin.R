##copy of 2016-09-16/rePlotCiberSort.R, adjusted for cibersort with dermals + skin

##replot cibersort
require(synapseClient)
synapseLogin()

cs.res<-read.delim(synGet('syn7810244')@filePath,sep="\t",header=T)

##the PBK ids are missing from table, so need to query annotations
res<-synQuery("select patientID,tissueID,sampleID from entity where parentId=='syn5492805'")
#map<-unique(res)

#from table get generic tumor id
tres<-synTableQuery("SELECT Patient,RnaID,TumorNumber,'RNASeq (Cufflinks)' FROM syn5556216 where RnaID is not NULL")@values

idx<-match(res$entity.id,tres$`RNASeq (Cufflinks)`)

dres<-res[which(!is.na(idx)),]
tres<-tres[idx[which(!is.na(idx))],]

full.map<-cbind(dres,tres)

source('../../bin/encodeSkinRNASeq.R')

sampleIds<-sapply(cs.res$Input.Sample,function(x){
  y=which(full.map$entity.sampleID==gsub('”','',x))
  paste("Patient",full.map$Patient[y],"Tumor",full.map$TumorNumber[y])
  
})
  
cibersort.df <- t(read.table(synGet("syn7810244")@filePath,sep='\t',header=T))
header <- c(cibersort.df[1,])
cibersort.df <- cibersort.df[2:26,]
colnames(cibersort.df) <- header
skin.type <- getSampleNamesForMatrix(cibersort.df)
skin.IDs <- rownames(skin.type)
skin.type.s <- c(paste((skin.type$Sample), skin.IDs, sep = "-"))

full.map$entity.sampleID <- gsub("-", ".", full.map$entity.sampleID, fixed = TRUE)
full.map$entity.sampleID <- paste("X",full.map$entity.sampleID, sep = "")

all.names <- c(sampleIds[1:33], skin.type.s[34:66])

#map tumors to sample ids

rownames(cs.res)<-all.names
cs.res<-cs.res[,-1]

this.file=''
library(pheatmap)
pheatmap(cs.res[,1:22],annotation_row=cs.res[,23:25],cellheight=10,cellwidth = 10)
synStore(File('ciberSortRePlotted.png',parentId='syn5809348'),executed=list(list(url=this.file)))
