source("../../bin/encodeSkinRNASeq.R")

## now get dermal NF data and cluster alongside
source("../../dermalNF/bin/dermalNFData.R")
dermals = rna_fpkm_matrix(byIsoform = FALSE)

# count_matrix(stored=TRUE,doNorm=FALSE,minCount=2,doLogNorm=FALSE,doVoomNorm=TRUE)
skin = getGeneNamesForMatrix(getEncodeSkinMatrix(metric = "FPKM", alignment = "hg19", doVoomNorm = FALSE))

over = intersect(rownames(dermals), rownames(skin))
## which annotation should we do? Are they really just duplicates of one another?

## step 1 - just combine all
comb = cbind(dermals[over, ], skin[over, ])

## step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr <- 1:nrow(comb)  #which(apply(comb,1,function(x) all(x>0.2)))
expr <- setdiff(expr, expr[union(grep("MIR", rownames(comb)[expr]), grep("SNO", rownames(comb)[expr]))])

## step 3, normalize
require(limma)
comb.norm = data.frame(voomWithQualityWeights(comb[expr, ])$E)
comb.norm$Gene = rownames(comb.norm)
comb.norm2 = comb.norm[1:66]


## write TSV for CIBERSORT - all samples, voom normalized
write.table(comb.norm2, file = "dNF_ENCODE.txt", sep = "\t", quote = FALSE, col.names = NA)

## write TSV for CIBERSORT - only skin (ENCODE) samples, before limma
write.table(skin, file = "ENCODE.txt", sep = "\t", quote = FALSE, col.names = NA)

