source("../../bin/encodeSkinRNASeq.R")

##now get dermal NF data and cluster alongisde
source("../../dermalNF/bin/dermalNFData.R")
dermals=rna_fpkm_matrix(byIsoform=FALSE)

##step 1 - just combine all
comb=dermals

##step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr<-1:nrow(comb) #which(apply(comb,1,function(x) all(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(comb)[expr]),grep("SNO",rownames(comb)[expr]))])

##step 3, normalize
require(limma)
comb.norm=data.frame(voomWithQualityWeights(comb[expr,])$E)
comb.norm$Gene=rownames(comb.norm)
comb.norm2=comb.norm


library(estimate)
write.table(comb.norm2, file = "combdata.txt", sep = "\t", quote = FALSE)
dNF.data <- c("combdata.txt")

filteredCommonGenes <- "combnorm.gct"
dNF.estimate <- "ESTIMATE.gct"

filterCommonGenes(input.f = dNF.data, output.f = filteredCommonGenes, id = "GeneSymbol")

estimateScore(filteredCommonGenes, dNF.estimate)

plotPurity(scores = dNF.estimate)

estimateScores.m <- t(dplyr::slice(read.delim(dNF.estimate)[2:35], 2:6))
estimateCols <- estimateScores.m[1:1,2:5]
estimateRows <- estimateScores.m[2:33,1:1]
estimateScores.m <- estimateScores.m[2:33, 2:5]

dims <- dim(estimateScores.m)
estimateScores.m <-as.numeric(estimateScores.m)
dim(estimateScores.m) <- dims

rownames(estimateScores.m) <- estimateRows
colnames(estimateScores.m) <- estimateCols

EstimateImmuneStromal <- estimateScores.m[1:32,1:3]
ImmuneStromal <- estimateScores.m[1:32,1:2]

library(pheatmap)
pheatmap(EstimateImmuneStromal, cellwidth = 8, cellheight = 8, cluster_cols = FALSE)
pheatmap(ImmuneStromal, cellwidth = 8, cellheight = 8, cluster_cols = FALSE)
