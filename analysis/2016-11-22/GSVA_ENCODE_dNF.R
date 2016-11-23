library(GSEABase)
library(GSVA)

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

## requires gene set collection in wd from MUTSIGDB
## http://software.broadinstitute.org/gsea/downloads.jsp
oncogenic.sigs <- getGmt("c6.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "c6"), 
    geneIdType = SymbolIdentifier())
hallmark.sigs <- getGmt("h.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"), 
    geneIdType = SymbolIdentifier())
immunologic.sigs <- getGmt("c7.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"), 
    geneIdType = SymbolIdentifier())

## oncogenic signature gene set collection GSVA
oncogenic.GSVA <- gsva(comb.norm2, oncogenic.sigs)$es.obs

## hallmarks signature gene set collection GSVA
hallmark.GSVA <- gsva(comb.norm2, hallmark.sigs)$es.obs

## immunologic signature gene set collection GSVA
immunologic.GSVA <- gsva(comb.norm2, immunologic.sigs)$es.obs

adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

library(pheatmap)
pheatmap(oncogenic.GSVA, fontsize_col = 9, fontsize_row = 1.3, border_color = FALSE)
pheatmap(hallmark.GSVA, fontsize_col = 9, fontsize_row = 4, border_color = FALSE)
pheatmap(immunologic.GSVA, fontsize_col = 9, fontsize_row = 0.25, border_color = FALSE)

