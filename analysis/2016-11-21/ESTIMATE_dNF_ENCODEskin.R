source("../../bin/encodeSkinRNASeq.R")

##now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals=rna_fpkm_matrix(byIsoform=FALSE)

#count_matrix(stored=TRUE,doNorm=FALSE,minCount=2,doLogNorm=FALSE,doVoomNorm=TRUE)
skin=getGeneNamesForMatrix(getEncodeSkinMatrix(metric='FPKM',alignment='hg19',doVoomNorm = FALSE))

over=intersect(rownames(dermals),rownames(skin))
##which annotation should we do? Are they really just duplicates of one another? 

##step 1 - just combine all
comb=cbind(dermals[over,],skin[over,])

##step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr<-1:nrow(comb) #which(apply(comb,1,function(x) all(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(comb)[expr]),grep("SNO",rownames(comb)[expr]))])

##step 3, normalize
require(limma)
comb.norm=data.frame(voomWithQualityWeights(comb[expr,])$E)
comb.norm$Gene=rownames(comb.norm)
comb.norm2=comb.norm[1:66]


library(estimate)
write.table(comb.norm2, file = "combdata.txt", sep = "\t", quote = FALSE)
dNF.skin <- c("combdata.txt")

filteredCommonGenes <- "combnorm.gct"
dNF.skin.estimate <- "ESTIMATE.gct"

filterCommonGenes(input.f = dNF.skin, output.f = filteredCommonGenes, id = "GeneSymbol")

estimateScore(filteredCommonGenes, dNF.skin.estimate)

plotPurity(scores = dNF.skin.estimate)

estimateScores.m <- t(dplyr::slice(read.delim(dNF.skin.estimate)[2:68], 2:6))
estimateCols <- estimateScores.m[1:1,2:5]
estimateRows <- estimateScores.m[2:66,1:1]
estimateScores.m <- estimateScores.m[2:66, 2:5]

dims <- dim(estimateScores.m)
estimateScores.m <-as.numeric(estimateScores.m)
dim(estimateScores.m) <- dims

rownames(estimateScores.m) <- estimateRows
colnames(estimateScores.m) <- estimateCols

EstimateImmuneStromal <- estimateScores.m[1:65,1:3]
ImmuneStromal <- estimateScores.m[1:65,1:2]

pheatmap:pheatmap(EstimateImmuneStromal, cellwidth = 8, cellheight = 8, cluster_cols = FALSE)
pheatmap:pheatmap(ImmuneStromal, cellwidth = 8, cellheight = 8, cluster_cols = FALSE)
