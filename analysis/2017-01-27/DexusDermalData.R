library(dexus)
library(pheatmap)

## now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals = rna_fpkm_matrix(byIsoform = FALSE)

comb = dermals

expr <- 1:nrow(comb)  ##which(apply(comb,1,function(x) all(x>0.2)))
expr <- setdiff(expr, expr[union(grep("MIR", rownames(comb)[expr]), grep("SNO", rownames(comb)[expr]))])

require(limma)
comb.norm = data.frame(voomWithQualityWeights(comb[expr, ])$E)
comb.norm$Gene = rownames(comb.norm)
comb.norm <- comb.norm[,1:33]
comb.norm <- as.matrix(comb.norm)

##the PBK ids are missing from table, so need to query annotations
res<-synQuery("select patientID,tissueID,sampleID from entity where parentId=='syn5492805'")
#map<-unique(res)
#from table get generic tumor id

tres<-synTableQuery("SELECT Patient,RnaID,TumorNumber,'RNASeq (Cufflinks)' FROM syn5556216 where RnaID is not NULL")@values
idx<-match(res$entity.id,tres$`RNASeq (Cufflinks)`)
dres<-res[which(!is.na(idx)),]
tres<-tres[idx[which(!is.na(idx))],]
full.map<-cbind(dres,tres)

full.map$entity.sampleID <- gsub("-", ".", full.map$entity.sampleID, fixed = TRUE)
full.map$entity.sampleID <- paste("X",full.map$entity.sampleID, sep = "")
sampleIds<-sapply(colnames(comb),function(x){
  y=which(full.map$entity.sampleID==gsub('”','',x))
  paste("patient",full.map$Patient[y],"tumor",full.map$TumorNumber[y])
})

colnames(comb.norm) <- sampleIds
comb.norm <- comb.norm[,1:33]

DermalDexus <- dexus(comb.norm[,1:33], normalization ='none')
informativetranscripts <- INI(DermalDexus, threshold=0.2)
plot(sort(informativetranscripts), idx = 1:30)
top30names<-(informativetranscripts@transcriptNames)[1:30]
top90names<-(informativetranscripts@transcriptNames)[1:90]

comb.norm2<-as.data.frame(comb.norm)
comb.norm2$Genes<-rownames(comb.norm2)

library(wesanderson)
wes<-wes_palette("Zissou", 200, type = "continuous")
  
heatmap<-filter(comb.norm2, Genes %in% top30names)
heatmap <- as.matrix(heatmap[,-34])
rownames(heatmap) <- top30names
pheatmap(heatmap, color = wes)

heatmap<-filter(comb.norm2, Genes %in% top90names)
heatmap <- as.matrix(heatmap[,-34])
rownames(heatmap) <- top90names
pheatmap(heatmap, fontsize_row = 5, color = wes)


