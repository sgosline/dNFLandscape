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

## oncogenic signature gene set collection GSVA and ssGSEA
oncogenic.GSVA <- gsva(comb.norm2, oncogenic.sigs)$es.obs
oncogenic.ssGSEA <- gsva(comb.norm2, oncogenic.sigs, method = "ssgsea")

## hallmarks signature gene set collection GSVA
hallmark.GSVA <- gsva(comb.norm2, hallmark.sigs)$es.obs
hallmark.ssGSEA <- gsva(comb.norm2, hallmark.sigs, method = "ssgsea")

## immunologic signature gene set collection GSVA
immunologic.GSVA <- gsva(comb.norm2, immunologic.sigs)$es.obs
immunologic.ssGSEA <- gsva(comb.norm2, immunologic.sigs, method = "ssgsea")

adjPvalueCutoff <- 0.001
logFCcutoff <- log2(2)

library(pheatmap)
pheatmap(oncogenic.GSVA, fontsize_col = 9, fontsize_row = 1.3, border_color = FALSE)
pheatmap(oncogenic.ssGSEA, fontsize_col = 9, fontsize_row = 1.3, border_color = FALSE)
pheatmap(hallmark.GSVA, fontsize_col = 9, fontsize_row = 4, border_color = FALSE)
pheatmap(hallmark.ssGSEA, fontsize_col = 9, fontsize_row = 4, border_color = FALSE)
pheatmap(immunologic.GSVA, fontsize_col = 9, fontsize_row = 0.25, border_color = FALSE)
pheatmap(immunologic.ssGSEA, fontsize_col = 9, fontsize_row = 4, border_color = FALSE)


##linear model to identify biggest oncogenic contributors to dNF (1) or skin (0)
oncogenic.ssGSEA.t<-t(oncogenic.ssGSEA)
oncogenic.ssGSEA.df <- as.data.frame(oncogenic.ssGSEA.t)

annotation <- c(rep.int(1, 33), rep.int(0, 33))

oncogenic1.lm2 <- lm(oncogenic.ssGSEA.df2$annotation ~ . , data = oncogenic.ssGSEA.df2[2:64])
summary(oncogenic1.lm2)
plot(oncogenic.lm2)

oncogenic2.lm2 <- lm(oncogenic.ssGSEA.df2$annotation ~ . , data = oncogenic.ssGSEA.df2[65:128])
summary(oncogenic2.lm2)
plot(oncogenic2.lm2)

oncogenic3.lm2 <- lm(oncogenic.ssGSEA.df2$annotation ~ . , data = oncogenic.ssGSEA.df2[129:190])
summary(oncogenic3.lm2)
plot(oncogenic3.lm2)

annotation.names <- c(rep.int("dNF", 33), rep.int("skin", 33))
oncogenic.ssGSEA.df3 <- cbind(annotation.names, oncogenic.ssGSEA.df)

plot(x= oncogenic.ssGSEA.df3$annotation.names, y = oncogenic.ssGSEA.df3$STK33_NOMO_DN)
t.test(oncogenic.ssGSEA.df3$STK33_NOMO_DN ~ oncogenic.ssGSEA.df3$annotation.names)

plot(x= oncogenic.ssGSEA.df3$annotation.names, y = oncogenic.ssGSEA.df3$GCNP_SHH_UP_LATE.V1_UP)
t.test(oncogenic.ssGSEA.df3$GCNP_SHH_UP_LATE.V1_UP ~ oncogenic.ssGSEA.df3$annotation.names)

plot(x= oncogenic.ssGSEA.df3$annotation.names, y = oncogenic.ssGSEA.df3$PIGF_UP.V1_DN)
t.test(oncogenic.ssGSEA.df3$PIGF_UP.V1_DN ~ oncogenic.ssGSEA.df3$annotation.names)

##

hallmark.GSVA.t<-t(hallmark.GSVA)
hallmark.GSVA.df <- as.data.frame(hallmark.GSVA.t)

hallmark.GSVA.df <- cbind(annotation, hallmark.GSVA.df)
hallmarkGSVA.lm <- lm(hallmark.GSVA.df$annotation ~ . , data = hallmark.GSVA.df[2:50])
summary(hallmarkGSVA.lm)
plot(hallmarkGSVA.lm)

hallmarkGSVA.lm2 <- lm(hallmark.GSVA.df$annotation ~ hallmark.GSVA.df$HALLMARK_HEDGEHOG_SIGNALING + hallmark.GSVA.df$HALLMARK_ESTROGEN_RESPONSE_LATE + 
                         hallmark.GSVA.df$HALLMARK_INTERFERON_ALPHA_RESPONSE + hallmark.GSVA.df$HALLMARK_P53_PATHWAY
                   + hallmark.GSVA.df$HALLMARK_APOPTOSIS + hallmark.GSVA.df$HALLMARK_KRAS_SIGNALING_UP, data = hallmark.GSVA.df)
summary(hallmarkGSVA.lm2)
plot(hallmarkGSVA.lm2)


##linear model to identify biggest hallmark contributors to dNF (1) or skin (0)

hallmark.ssGSEA.t<-t(hallmark.ssGSEA)
hallmark.ssGSEA.df <- as.data.frame(hallmark.ssGSEA.t)

hallmark.ssGSEA.df2 <- cbind(annotation, hallmark.ssGSEA.df)
hallmark.lm <- lm(hallmark.ssGSEA.df2$annotation ~ . , data = hallmark.ssGSEA.df2)
summary(hallmark.lm)
plot(hallmark.lm)

annotation.names <- c(rep.int("dNF", 33), rep.int("skin", 33))
hallmark.ssGSEA.df3 <- cbind(annotation.names, hallmark.ssGSEA.df)

plot(x= hallmark.ssGSEA.df3$annotation.names, y = hallmark.ssGSEA.df3$HALLMARK_HEDGEHOG_SIGNALING)
t.test(hallmark.ssGSEA.df3$HALLMARK_HEDGEHOG_SIGNALING ~ hallmark.ssGSEA.df3$annotation.names)

plot(x= hallmark.ssGSEA.df3$annotation.names, y = hallmark.ssGSEA.df3$HALLMARK_ESTROGEN_RESPONSE_LATE)
t.test(hallmark.ssGSEA.df3$HALLMARK_ESTROGEN_RESPONSE_LATE ~ hallmark.ssGSEA.df3$annotation.names)

plot(x= hallmark.ssGSEA.df3$annotation.names, y = hallmark.ssGSEA.df3$HALLMARK_INTERFERON_ALPHA_RESPONSE)
t.test(hallmark.ssGSEA.df3$HALLMARK_INTERFERON_ALPHA_RESPONSE ~ hallmark.ssGSEA.df3$annotation.names)

plot(x= hallmark.ssGSEA.df3$annotation.names, y = hallmark.ssGSEA.df3$HALLMARK_P53_PATHWAY)
t.test(hallmark.ssGSEA.df3$HALLMARK_P53_PATHWAY ~ hallmark.ssGSEA.df3$annotation.names)

hallmark.lm2 <- lm(hallmark.ssGSEA.df2$annotation ~ hallmark.ssGSEA.df2$HALLMARK_HEDGEHOG_SIGNALING + hallmark.ssGSEA.df2$HALLMARK_ESTROGEN_RESPONSE_LATE + 
   hallmark.ssGSEA.df2$HALLMARK_INTERFERON_ALPHA_RESPONSE + hallmark.ssGSEA.df2$HALLMARK_P53_PATHWAY
   + hallmark.ssGSEA.df2$HALLMARK_APOPTOSIS + hallmark.ssGSEA.df2$HALLMARK_KRAS_SIGNALING_UP, data = hallmark.ssGSEA.df2)
summary(hallmark.lm2)
plot(hallmark.lm2)

