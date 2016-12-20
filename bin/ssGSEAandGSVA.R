library(GSEABase)
library(GSVA)
library(ggplot2)
source("../../bin/encodeSkinRNASeq.R")

## now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals = rna_fpkm_matrix(byIsoform = FALSE)

# count_matrix(stored=TRUE,doNorm=FALSE,minCount=2,doLogNorm=FALSE,doVoomNorm=TRUE)
skin = getGeneNamesForMatrix(getEncodeSkinMatrix(metric = "FPKM", alignment = "hg19",
                                                 doVoomNorm = FALSE))
over = intersect(rownames(dermals), rownames(skin))

## which annotation should we do? Are they really just duplicates of one another?

## step 1 - just combine all
comb = cbind(dermals[over, ], skin[over, ])

## step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr <- 1:nrow(comb)  ##which(apply(comb,1,function(x) all(x>0.2)))
expr <- setdiff(expr, expr[union(grep("MIR", rownames(comb)[expr]), grep("SNO", rownames(comb)[expr]))])

## step 3, normalize
require(limma)
comb.norm = data.frame(voomWithQualityWeights(comb[expr, ])$E)
comb.norm$Gene = rownames(comb.norm)
comb.norm2 = comb.norm[1:66]
comb.norm2 <- as.matrix(comb.norm2)
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
cibersort.df <- t(read.table(synGet("syn7810244")@filePath,sep='\t',header=T))
header <- c(cibersort.df[1,])
cibersort.df <- cibersort.df[2:26,]
colnames(cibersort.df) <- header
skin.type <- getSampleNamesForMatrix(cibersort.df)
skin.IDs <- rownames(skin.type)
skin.type.s <- c(paste((skin.type$Sample), skin.IDs, sep = " "))
full.map$entity.sampleID <- gsub("-", ".", full.map$entity.sampleID, fixed = TRUE)
full.map$entity.sampleID <- paste("X",full.map$entity.sampleID, sep = "")
sampleIds<-sapply(cs.res$Input.Sample,function(x){
  y=which(full.map$entity.sampleID==gsub('”','',x))
  paste("patient",full.map$Patient[y],"tumor",full.map$TumorNumber[y])
})
all.names <- c(sampleIds[1:33], skin.type.s[34:66])
colnames(comb.norm2) <- all.names

## requires gene set collection in wd from MUTSIGDB
## http://software.broadinstitute.org/gsea/downloads.jsp
#oncogenic.sigs <- getGmt("../../data/c6.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "c6"),
#                         geneIdType = SymbolIdentifier())
hallmark.sigs <- getGmt("../../data/h.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"),
                        geneIdType = SymbolIdentifier())
#immunologic.sigs <- getGmt("../../data/c7.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"),
#                          geneIdType = SymbolIdentifier())

## oncogenic signature gene set collection sGSEA
#oncogenic.ssGSEA <- gsva(comb.norm2, oncogenic.sigs, method = "ssgsea")
## hallmarks signature gene set collection ssgsea
hallmark.ssGSEA <- gsva(comb.norm2, hallmark.sigs, method = "ssgsea")
## immunologic signature gene set collection ssgsea
#immunologic.ssGSEA <- gsva(comb.norm2, immunologic.sigs, method = "ssgsea")