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
oncogenic.sigs <- getGmt("../../data/c6.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "c6"),
                         geneIdType = SymbolIdentifier())
hallmark.sigs <- getGmt("../../data/h.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"),
                        geneIdType = SymbolIdentifier())
#immunologic.sigs <- getGmt("../../data/c7.all.v5.2.symbols.gmt", collectionType = BroadCollection(category = "h"),
#                           geneIdType = SymbolIdentifier())

## oncogenic signature gene set collection sGSEA
oncogenic.ssGSEA <- gsva(comb.norm2, oncogenic.sigs, method = "ssgsea")
## hallmarks signature gene set collection ssgsea
hallmark.ssGSEA <- gsva(comb.norm2, hallmark.sigs, method = "ssgsea")
## immunologic signature gene set collection ssgsea
immunologic.ssGSEA <- gsva(comb.norm2, immunologic.sigs, method = "ssgsea")

oncogenic.ssGSEA.t<-t(oncogenic.ssGSEA)
annotation <- as.factor(c(rep.int(1, 33), rep.int(0, 33)))

oncogenic.ssGSEA.df <- as.data.frame(oncogenic.ssGSEA.t)
oncogenic.ssGSEA.df <- cbind(annotation, oncogenic.ssGSEA.df)

##modified from http://dni-institute.in/blogs/random-forest-using-r-step-by-step-tutorial/
sample.ind <- sample(2, 
                     nrow(oncogenic.ssGSEA.df),
                     replace = T,
                     prob = c(0.6,0.4))
oncogenic.ssGSEA.dev <- oncogenic.ssGSEA.df[sample.ind==1,]
oncogenic.ssGSEA.val <- oncogenic.ssGSEA.df[sample.ind==2,]

library(randomForest)
oncogenic.rF<-randomForest(formula = oncogenic.ssGSEA.dev$annotation ~ ., data = oncogenic.ssGSEA.dev)

summary(oncogenic.rF)
plot(oncogenic.rF)
varImpPlot(oncogenic.rF, n.var=20)

predict(oncogenic.rF, oncogenic.ssGSEA.dev)
predict(oncogenic.rF, oncogenic.ssGSEA.val)

oncogenic.ssGSEA.t<-t(oncogenic.ssGSEA)
annotation <- as.factor(c(rep.int(1, 33), rep.int(0, 33)))

oncogenic.ssGSEA.df <- as.data.frame(oncogenic.ssGSEA.t)
oncogenic.ssGSEA.df <- cbind(annotation, oncogenic.ssGSEA.df)

##modified from http://dni-institute.in/blogs/random-forest-using-r-step-by-step-tutorial/
#sample.ind <- sample(2, 
#                     nrow(oncogenic.ssGSEA.df),
#                     replace = T,
#                     prob = c(0.6,0.4))
#oncogenic.ssGSEA.dev <- oncogenic.ssGSEA.df[sample.ind==1,]
#oncogenic.ssGSEA.val <- oncogenic.ssGSEA.df[sample.ind==2,]

oncogenic.rFull<-randomForest(formula = oncogenic.ssGSEA.df$annotation ~ ., data = oncogenic.ssGSEA.df)

sink("oncogenic.rFull.txt", append=FALSE)
summary(oncogenic.rFull)
sink()

pdf("oncogenic.randomforest.pdf")
plot(oncogenic.rFull)
treesize(oncogenic.rFull)
varImpPlot(oncogenic.rFull, n.var=20)
dev.off()

library(glmnet)
oncogenic.input <- oncogenic.ssGSEA.t
annotation <- c(rep.int(1, 33), rep.int(0, 33))


oncogenic.glm<-glmnet(x = oncogenic.input, y = annotation, alpha =1)
oncogenic.cvglm<-cv.glmnet(x = oncogenic.input, y = annotation, alpha =1)

plot.cv.glmnet(oncogenic.cvglm)
coef.cv.glmnet(oncogenic.cvglm)

pdf("oncogenicGLMnet.pdf")
plot.glmnet(oncogenic.glm)
dev.off()

sink("oncogenic.glm.txt", append=FALSE)
print.glmnet(oncogenic.glm)
coef.glmnet(oncogenic.glm)
sink()

#####################################################################
hallmark.ssGSEA.t<-t(hallmark.ssGSEA)
annotation <- as.factor(c(rep.int(1, 33), rep.int(0, 33)))

hallmark.ssGSEA.df <- as.data.frame(hallmark.ssGSEA.t)
hallmark.ssGSEA.df <- cbind(annotation, hallmark.ssGSEA.df)

##modified from http://dni-institute.in/blogs/random-forest-using-r-step-by-step-tutorial/
sample.ind <- sample(2, 
                    nrow(hallmark.ssGSEA.df),
                    replace = T,
                    prob = c(0.6,0.4))
hallmark.ssGSEA.dev <- hallmark.ssGSEA.df[sample.ind==1,]
annotations <- as.vector(hallmark.ssGSEA.dev[1], mode='numeric')
hallmark.ssGSEA.dev <- hallmark.ssGSEA.dev[-1]
hallmark.ssGSEA.val <- hallmark.ssGSEA.df[sample.ind==2,]
hallmark.ssGSEA.val <- hallmark.ssGSEA.val[1]

hall.dev <-cv.glmnet(x = hallmark.ssGSEA.dev, y = annotations, alpha =1)
predict.cv.glmnet(hall.dev, hallmark.ssGSEA.val, s="lamda.lse")



hallmark.rF<-randomForest(formula = hallmark.ssGSEA.df$annotation ~ ., data = hallmark.ssGSEA.df)

summary(hallmark.rF)

pdf("hallmark.randomforest.pdf")
plot(hallmark.rF)
varImpPlot(hallmark.rF, n.var=20)
dev.off()

sink("hallmark.rF.txt", append=FALSE)
summary(hallmark.rF)
sink()

#predict(hallmark.rF, hallmark.ssGSEA.dev)
#predict(hallmark.rF, hallmark.ssGSEA.val)

hallmark.input <- hallmark.ssGSEA.t
annotation <- c(rep.int(1, 33), rep.int(0, 33))

hallmark.glm<-glmnet(x = hallmark.input, y = annotation, alpha =1)
hallmark.cvglm<-cv.glmnet(x = hallmark.input, y = annotation, alpha =1)

hallmark.cvglm<-cv.glmnet(x = hallmark.input, y = annotation, alpha =1)


plot.cv.glmnet(hallmark.cvglm)
coef.cv.glmnet(hallmark.cvglm)
coef.glmnet(hallmark.glm)

pdf("hallmarkGLMnet.pdf")
plot.glmnet(hallmark.glm)
dev.off()

sink("hallmark.glm.txt", append=FALSE)
print.glmnet(hallmark.glm)
coef.glmnet(hallmark.glm)
sink()

#####################################################################

immunologic.ssGSEA.t<-t(immunologic.ssGSEA)
annotation <- as.factor(c(rep.int(1, 33), rep.int(0, 33)))

immunologic.ssGSEA.df <- as.data.frame(immunologic.ssGSEA.t)
immunologic.ssGSEA.df <- cbind(annotation, immunologic.ssGSEA.df)

##modified from http://dni-institute.in/blogs/random-forest-using-r-step-by-step-tutorial/
#sample.ind <- sample(2, 
#                     nrow(immunologic.ssGSEA.df),
        #             replace = T,
#                     prob = c(0.6,0.4))
#immunologic.ssGSEA.dev <- immunologic.ssGSEA.df[sample.ind==1,]
#immunologic.ssGSEA.val <- immunologic.ssGSEA.df[sample.ind==2,]

immunologic.ssGSEA.df2 <- immunologic.ssGSEA.df[1:411]

library(randomForest)
immunologic.rF<-randomForest(formula = immunologic.ssGSEA.df2$annotation ~ ., data = immunologic.ssGSEA.df2)

pdf("immunologic.randomforest.pdf")
plot(immunologic.rF)
varImpPlot(immunologic.rF, n.var=20)
dev.off()

sink("immunologic.rF.txt", append=FALSE)
summary(immunologic.rF)
sink()

#predict(immunologic.rF, immunologic.ssGSEA.dev)
#predict(immunologic.rF, immunologic.ssGSEA.val)

immunologic.input <- immunologic.ssGSEA.t
annotation <- c(rep.int(1, 33), rep.int(0, 33))

immunologic.glm<-glmnet(x = immunologic.input, y = annotation, alpha =1)
pdf("immunologicGLMnet.pdf")
plot.glmnet(immunologic.glm)
dev.off()

sink("immunologic.glm.txt", append=FALSE)
print.glmnet(immunologic.glm)
coef.glmnet(immunologic.glm)
sink()


#####################################################################
annotation.names <- c(rep.int("dNF", 33), rep.int("skin", 33))
oncogenic.ssGSEA.df2 <- cbind(annotation.names, oncogenic.ssGSEA.df[,-1])


#######################
hallmark.ssGSEA.df2 <- cbind(annotation.names, hallmark.ssGSEA.df[,-1])

hallmark.ttest <- as.data.frame(sapply(hallmark.ssGSEA.df2[2:51],
            function(x) t.test(x ~ hallmark.ssGSEA.df2$annotation.names)))

p <- c(unlist(hallmark.ttest[3,]))
p <- signif(p.adjust(p, method = "BH"), 3)
p.df <- as.data.frame(p)
rownames(p.df) <- sub(".p.value", "", rownames(p.df))
p.df <- p.df[order(p.df$p), , drop=FALSE]

write.table(p.df, file = 'hallmark_p_values.txt', sep = "\t")

pdf('hallmark_res.pdf')
apply(as.array(colnames(hallmark.ssGSEA.df2)[2:51]),1,FUN=function (x)
  qplot(annotation.names, hallmark.ssGSEA.df2[,x], data = hallmark.ssGSEA.df2) +
    geom_boxplot(fill=c("red", "blue")) +
    xlab("sample type") +
    ylab(x) +
    ggtitle(paste(x))
)
dev.off()

p.df$pathway <- rownames(p.df)
psig <- dplyr::filter(p.df,p<0.05)


######################################
oncogenic.ssGSEA.df2 <- cbind(annotation.names, oncogenic.ssGSEA.df[,-1])

oncogenic.ttest <- as.data.frame(sapply(oncogenic.ssGSEA.df2[2:190],
                                       function(x) t.test(x ~ oncogenic.ssGSEA.df2$annotation.names)))

p <- c(unlist(oncogenic.ttest[3,]))
p <- signif(p.adjust(p, method = "BH"), 3)
p.df <- as.data.frame(p)
rownames(p.df) <- sub(".p.value", "", rownames(p.df))
p.df <- p.df[order(p.df$p), , drop=FALSE]

write.table(p.df, file = 'oncogenic_p_values.txt', sep = "\t")




pdf('oncogenic_res.pdf')
apply(as.array(colnames(oncogenic.ssGSEA.df2)[2:190]),1,FUN=function (x)
  qplot(annotation.names, oncogenic.ssGSEA.df2[,x], data = oncogenic.ssGSEA.df2) +
    geom_boxplot(fill=c("red", "blue")) +
    xlab("sample type") +
    ylab(x) +
    ggtitle(paste(x))
)
dev.off()

p.df$pathway <- rownames(p.df)
psig <- dplyr::filter(p.df,p<0.05)

##############################################
#immunologic.ssGSEA.df2 <- cbind(annotation.names, immunologic.ssGSEA.df[,-1])

#immunologic.ttest <- as.data.frame(sapply(immunologic.ssGSEA.df2[2:51],
#                                        function(x) t.test(x ~ immunologic.ssGSEA.df2$annotation.names)))

#p <- c(unlist(immunologic.ttest[3,]))
#p <- signif(p.adjust(p, method = "BH"), 3)
#p.df <- as.data.frame(p)
#rownames(p.df) <- sub(".p.value", "", rownames(p.df))

#write.table(p.df, file = 'immunologic_p_values.txt', sep = "\t")

#pdf('immuno_res.pdf')
#apply(as.array(colnames(immunologic.ssGSEA.df2)[2:4873]),1,FUN=function (x)
#  qplot(annotation.names, immunologic.ssGSEA.df2[,x], data = immunologic.ssGSEA.df2) +
#    geom_boxplot(fill=c("red", "blue")) +
#    xlab("sample type") +
#    ylab("ssGSEA output") +
#    ggtitle(paste(x))
#)
#dev.off()

#p.df$pathway <- rownames(p.df)
#psig <- dplyr::filter(p.df,p<0.05)