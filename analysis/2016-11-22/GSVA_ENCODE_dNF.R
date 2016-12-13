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
  paste("Patient",full.map$Patient[y],"Tumor",full.map$TumorNumber[y])
})
all.names <- c(sampleIds[1:33], skin.type.s[34:66])
colnames(comb.norm2) <- all.names

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
this.file <- "https://github.com/allaway/dNFLandscape/blob/master/analysis/2016-11-22/GSVA_ENCODE_dNF.R"
pheatmap(oncogenic.GSVA, fontsize_col = 2, fontsize_row = 1.3, border_color = FALSE, filename = 'oncogenic.GSVA.png')
synStore(File('oncogenic.GSVA.png', parentId = 'syn7818711'), used = 'syn7818712', executed = this.file, activityName = 'GSVA')
pheatmap(oncogenic.ssGSEA, fontsize_col = 2, fontsize_row = 1.3, border_color = FALSE, filename = 'oncogenic.ssGSEA.png')
synStore(File('oncogenic.ssGSEA.png', parentId = 'syn7818711'), used = 'syn7818712', executed = this.file, activityName = 'ssGSEA')
pheatmap(hallmark.GSVA, fontsize_col = 2, fontsize_row = 4, border_color = FALSE, filename = 'hallmark.GSVA.png')
synStore(File('hallmark.GSVA.png', parentId = 'syn7818711'), used = 'syn7818714', executed = this.file, activityName = 'GSVA')
pheatmap(hallmark.ssGSEA, fontsize_col = 2, fontsize_row = 4, border_color = FALSE, filename = 'hallmark.ssGSEA.png')
synStore(File('hallmark.ssGSEA.png', parentId = 'syn7818711'), used = 'syn7818714', executed = this.file, activityName = 'ssGSEA')
pheatmap(immunologic.GSVA, fontsize_col = 2, fontsize_row = 0.25, border_color = FALSE, filename = 'immunologic.GSVA.png')
synStore(File('immunologic.GSVA.png', parentId = 'syn7818711'), used = 'syn7818713', executed = this.file, activityName = 'GSVA')
pheatmap(immunologic.ssGSEA, fontsize_col = 2, fontsize_row = 0.25, border_color = FALSE, filename = 'immunologic.ssGSEA.png')
synStore(File('immunologic.ssGSEA.png', parentId = 'syn7818711'), used = 'syn7818713', executed = this.file, activityName = 'ssGSEA')

##linear model to identify biggest oncogenic contributors to dNF (1) or skin (0)
oncogenic.ssGSEA.t<-t(oncogenic.ssGSEA)
oncogenic.ssGSEA.df <- as.data.frame(oncogenic.ssGSEA.t)
annotation <- c(rep.int(1, 33), rep.int(0, 33))
oncogenic.ssGSEA.df <- cbind(oncogenic.ssGSEA.df, annotation)
oncogenic1.lm2 <- lm(oncogenic.ssGSEA.df$annotation ~ . , data = oncogenic.ssGSEA.df[2:64])
sink('oncogenic1.lm.txt')
summary(oncogenic1.lm2, file = "oncogenic1.lm.txt")
sink()
synStore(File('oncogenic1.lm.txt', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')
png("oncogenic1.lm.png")
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(oncogenic1.lm2, ask = FALSE)
dev.off()
synStore(File('oncogenic1.lm.png', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')
oncogenic2.lm2 <- lm(oncogenic.ssGSEA.df$annotation ~ . , data = oncogenic.ssGSEA.df[65:128])
sink('oncogenic2.lm.txt')
summary(oncogenic2.lm2, file= "oncogenic2.lm.txt")
sink()
synStore(File('oncogenic2.lm.txt', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')
png("oncogenic2.lm.png")
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(oncogenic2.lm2, ask = FALSE)
dev.off()
synStore(File('oncogenic2.lm.png', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')

oncogenic3.lm2 <- lm(oncogenic.ssGSEA.df$annotation ~ . , data = oncogenic.ssGSEA.df[129:190])
sink('oncogenic3.lm.txt')
summary(oncogenic3.lm2)
sink()
synStore(File('oncogenic3.lm.txt', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')

png("oncogenic3.lm.png")
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(oncogenic3.lm2, ask = FALSE)
dev.off()
synStore(File('oncogenic3.lm.png', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')
annotation.names <- c(rep.int("dNF", 33), rep.int("skin", 33))
oncogenic.ssGSEA.df2 <- cbind(annotation.names, oncogenic.ssGSEA.df)


pdf('OncogenicBoxPlots.pdf')

##plot significant contributors to lm(oncogenic) as individual boxplots + ttest to determine if different
p<-signif(t.test(oncogenic.ssGSEA.df2$EGFR_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, EGFR_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
p<-signif(t.test(oncogenic.ssGSEA.df2$GCNP_SHH_UP_EARLY.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, GCNP_SHH_UP_EARLY.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$GCNP_SHH_UP_LATE.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, GCNP_SHH_UP_LATE.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$CSR_LATE_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, CSR_LATE_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))

##PIGF is a non-essential angiogenesis factor that has been proposed as a therapeutic target
##for medulloblastoma among other cancers - appears to be differentially regulated in the two sample sets
p<-signif(t.test(oncogenic.ssGSEA.df2$PIGF_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, PIGF_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$P53_DN.V2_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, P53_DN.V2_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$E2F3_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, E2F3_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$SRC_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, SRC_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$CAHOY_OLIGODENDROCUTIC~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, CAHOY_OLIGODENDROCUTIC, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$RPS14_DN.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, RPS14_DN.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$IL15_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, IL15_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$IL15_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, IL15_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$PDGF_ERK_DN.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, PDGF_ERK_DN.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$PDGF_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, PDGF_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$YAP1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, YAP1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$YAP1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, YAP1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$SIRNA_EIF4GI_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, SIRNA_EIF4GI_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$STK33_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, STK33_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$STK33_NOMO_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, STK33_NOMO_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$STK33_SKM_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, STK33_SKM_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.DF.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.DF.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$TBK1.DN.48HRS_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, TBK1.DN.48HRS_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.300_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.300_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.600_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.600_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.600.LUNG.BREAST_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.600.LUNG.BREAST_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.BREAST_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.BREAST_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.BREAST_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.BREAST_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.KIDNEY_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.KIDNEY_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.LUNG_UP.V1_DN~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.LUNG_UP.V1_DN, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.LUNG_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.LUNG_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.LUNG.BREAST_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.LUNG.BREAST_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$KRAS.PROSTATE_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, KRAS.PROSTATE_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
p<-signif(t.test(oncogenic.ssGSEA.df2$LEF1_UP.V1_UP~ oncogenic.ssGSEA.df2$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 33)
qplot(annotation.names, LEF1_UP.V1_UP, data = oncogenic.ssGSEA.df2) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =",p))
dev.off()

synStore(File("OncogenicBoxPlots.pdf", parentId = 'syn7818711'), used = 'syn7818712', executed = this.file)

hallmark.GSVA.t<-t(hallmark.GSVA)
hallmark.GSVA.df <- as.data.frame(hallmark.GSVA.t)
hallmark.GSVA.df <- cbind(annotation, hallmark.GSVA.df)
hallmarkGSVA.lm <- lm(hallmark.GSVA.df$annotation ~ . , data = hallmark.GSVA.df[2:50])
summary(hallmarkGSVA.lm)
sink("hallmarkGSVA.lm.txt")
summary(hallmarkGSVA.lm)
sink()
synStore(File('hallmarkGSVA.lm.txt', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')

png("hallmarkGSVA.lm.png")
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(hallmarkGSVA.lm, ask = FALSE)
dev.off()
synStore(File('hallmarkGSVA.lm.png', parentId = 'syn7818711'),
         used = 'syn7818714',
         executed = this.file,
         activityName = 'linear modeling')

##linear model to identify biggest hallmark contributors to dNF (1) or skin (0)
hallmark.ssGSEA.t<-t(hallmark.ssGSEA)
hallmark.ssGSEA.df <- as.data.frame(hallmark.ssGSEA.t)
hallmark.ssGSEA.df2 <- cbind(annotation, hallmark.ssGSEA.df)
hallmark.lm <- lm(hallmark.ssGSEA.df2$annotation ~ . , data = hallmark.ssGSEA.df2)

sink("hallmarkssGSEA.lm.txt")
summary(hallmark.lm)
sink()
synStore(File('hallmarkssGSEA.lm.txt', parentId = 'syn7818711'),
         used = 'syn7818712',
         executed = this.file,
         activityName = 'linear modeling')

png("hallmarkssGSEA.lm.png")
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
plot(hallmark.lm, ask = FALSE)
dev.off()
synStore(File('hallmarkssGSEA.lm.png', parentId = 'syn7818711'),
         used = 'syn7818714',
         executed = this.file,
         activityName = 'linear modeling')
annotation.names <- c(rep.int("dNF", 33), rep.int("skin", 33))
hallmark.ssGSEA.df3 <- cbind(annotation.names, hallmark.ssGSEA.df)
pdf("HallmarkBoxPlots.pdf")
p<-signif(t.test(hallmark.ssGSEA.df3$HALLMARK_HEDGEHOG_SIGNALING ~ hallmark.ssGSEA.df3$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 4)
qplot(annotation.names, HALLMARK_HEDGEHOG_SIGNALING, data = hallmark.ssGSEA.df3) +
  geom_boxplot(fill=c("blue", "red"))  +
  ggtitle(paste("FDR adjusted p-value =", p))
p<-signif(t.test(hallmark.ssGSEA.df3$HALLMARK_ESTROGEN_RESPONSE_LATE ~ hallmark.ssGSEA.df3$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 4)
qplot(annotation.names, HALLMARK_ESTROGEN_RESPONSE_LATE, data = hallmark.ssGSEA.df3) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =", p))
p<-signif(t.test(hallmark.ssGSEA.df3$HALLMARK_INTERFERON_ALPHA_RESPONSE ~ hallmark.ssGSEA.df3$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 4)
qplot(annotation.names, HALLMARK_INTERFERON_ALPHA_RESPONSE, data = hallmark.ssGSEA.df3) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =", p))
p<-signif(t.test(hallmark.ssGSEA.df3$HALLMARK_P53_PATHWAY ~ hallmark.ssGSEA.df3$annotation.names)$p.value, digits = 3)
p<-p.adjust(p, method = "BH", 4)
qplot(annotation.names, HALLMARK_P53_PATHWAY, data = hallmark.ssGSEA.df3) +
  geom_boxplot(fill=c("blue", "red")) +
  ggtitle(paste("FDR adjusted p-value =", p))
dev.off()
synStore(File("HallmarkBoxPlots.pdf", parentId = 'syn7818711'), used = 'syn7818714', executed = this.file)







