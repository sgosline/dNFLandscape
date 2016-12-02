source("../../bin/encodeSkinRNASeq.R")

##first check sanity of data from encode

#for(t in c('TPM','FPKM'))
#  for(a in c('hg19','grch38'))
#    clusterSamples(t,a)


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

require(tidyr)
norm.df<-tidyr::gather(comb.norm,"Sample","Expression",1:(ncol(comb.norm)-1))

mastsig <- c("KIT", "IL1RL1", "FCER1A", "MS4A2", "ENPP3", "HDC", "TPSAB1", "TPSB2", "TPSD1", "CMA1", "CPA3", "CTSG", "HPGD5", "LTC45")
#established mast cell expression signature from http://www.bloodjournal.org/content/123/17/e58.long
mastsig.df <- dplyr::filter(norm.df, Gene %in% mastsig)
mastsig.df <- tidyr::spread(mastsig.df, Gene, Expression)

mastsig.col <- mastsig.df[,-1]
rownames(mastsig.col) <- mastsig.df[,1]
mastsig.m <- as.matrix(mastsig.col)

samples <- rownames(mastsig.m)
sizes <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,13,14,5,7,11,16,20,15,8,5,10,7,8,6,15,10,15,13,10,10,13,10,12,5,11,20,10,3,18,12,20,13)
tumorsize.m <- matrix(data = sizes, nrow = 66, ncol = 1, dimnames = (list(samples, "Size")))
tumorsize.df <- as.data.frame(tumorsize.m)

pheatmap(mastsig.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)

mastsig2 <- c("MRGPRX2", "RGS13", "C1orf150", "SIGLEC6", "ERVFRD-1", "SVOPL", "C20orf118", "VWASA")
#new mast cell expression signature from http://www.bloodjournal.org/content/123/17/e58.long, many of these dont appear to be in this data set - perhaps different identfiers?

mastsig2.df <- dplyr::filter(norm.df, Gene %in% mastsig2)
mastsig2.df <- tidyr::spread(mastsig2.df, Gene, Expression)

mastsig2.col <- mastsig2.df[,-1]
rownames(mastsig2.col) <- mastsig2.df[,1]
mastsig2.m <- as.matrix(mastsig2.col)

pheatmap(mastsig2.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)

##linear model to determine if signatures are good predictors 

mastsig.df2 <- as.data.frame(cbind(sizes, mastsig.m))
mastsig.df2 <- dplyr::slice(mastsig.df2, 34:66)
rownames(mastsig.df2) <- samples[34:66]
mastsig.lm <- lm(mastsig.df2$sizes ~ ., data = mastsig.df2)
summary(mastsig.lm)

mastsig2.df2 <- as.data.frame(cbind(sizes, mastsig2.m))
mastsig2.df2 <- dplyr::slice(mastsig2.df2, 34:66)
rownames(mastsig2.df2) <- samples[34:66]
mastsig2.lm <- lm(mastsig2.df2$sizes ~ ., data = mastsig2.df2)
summary(mastsig2.lm)

plot(y = mastsig.df2$CMA1, x= mastsig.df2$sizes)

##handpicked candidates for lm()
mastsig3.df <- as.data.frame(cbind(sizes, CMA1=mastsig.df$CMA1, ENPP3=mastsig.df$ENPP3, `ERVFRD-1`=mastsig2.df$`ERVFRD-1`, 
  MRGPRX2=mastsig2.df$MRGPRX2))

mastsig3.df <- dplyr::slice(mastsig3.df, 34:66)
rownames(mastsig3.df) <- samples[34:66]
mastsig3.lm <- lm(mastsig3.df$sizes ~ . , data = mastsig3.df)
summary(mastsig3.lm)
plot(mastsig3.lm)

