library(GSEABase)
library(GSVA)
library(ggplot2)
library(dplyr)
library(viridis)

source("../../dermalNF/bin/dermalNFData.R")
dermals = rna_fpkm_matrix(byIsoform = FALSE)

## step 1 - just combine all
comb = dermals

## step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr <- 1:nrow(comb)  ##which(apply(comb,1,function(x) all(x>0.2)))
expr <- setdiff(expr, expr[union(grep("MIR", rownames(comb)[expr]), grep("SNO", rownames(comb)[expr]))])

## step 3, normalize
require(limma)
comb.norm = data.frame(voomWithQualityWeights(comb[expr, ])$E)
comb.norm$Gene = rownames(comb.norm)
comb.norm2 = comb.norm[1:33]
comb.norm2 <- as.matrix(comb.norm2)

colnames(comb.norm2) <- gsub(".", "-", colnames(comb.norm2), fixed = TRUE)
colnames(comb.norm2) <- gsub("X", "", colnames(comb.norm2), fixed = TRUE)

rna.id<-read.table(synGet("syn6156140")@filePath, sep = ",", header = TRUE)
pt.id<-synTableQuery("SELECT 'RNASeq (Cufflinks)', sampleIdentifier FROM syn5556216")@values %>% filter(!is.na(`RNASeq (Cufflinks)`))
colnames(pt.id) <- c("SynapseID", "patientID")
ids <- full_join(rna.id,pt.id)

colnames(comb.norm2)<-ids$patientID[colnames(comb.norm2)==ids$Sample.Id]
## requires gene set collection in wd from MUTSIGDB
## http://software.broadinstitute.org/gsea/downloads.jsp
hallmark.sigs <- getGmt("h.all.v6.0.symbols.gmt", collectionType = BroadCollection(category = "h"),
                        geneIdType = SymbolIdentifier())

this.file <- "https://github.com/allaway/dNFLandscape/blob/master/analysis/2016-11-22/GSVA_ENCODE_dNF_mini.R"

## hallmarks signature gene set collection GSVA
hallmark.ssGSEA <- as.data.frame(gsva(comb.norm2, hallmark.sigs, rnaseq=TRUE, method = "ssgsea"))
library(tidyr)
library(grid)
hallmark.ssGSEA$pathway<- rownames(hallmark.ssGSEA)
hallmark.ssGSEA <- tidyr::gather(hallmark.ssGSEA, patient, measurement, 1:33)
hallmark.ssGSEA$pathway <- sub("HALLMARK_", "", hallmark.ssGSEA$pathway)
hallmark.ssGSEA$pathway <- sub("_", " ", hallmark.ssGSEA$pathway)
hallmark.ssGSEA$pathway <- factor(hallmark.ssGSEA$pathway, levels = levels(reorder(hallmark.ssGSEA$pathway, -hallmark.ssGSEA$measurement, FUN = median)))

hallmark.ssGSEA$mean <- ave(hallmark.ssGSEA$measurement, as.factor(hallmark.ssGSEA$pathway), FUN=mean)

ggplot(hallmark.ssGSEA, aes(x=pathway, y=measurement, fill = mean)) +
  geom_boxplot() + 
  scale_fill_viridis(option = "B") +
  theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
  theme(plot.margin=unit(c(1,1,1,1), "cm")) +
  labs(x = "Hallmark Gene Set", y = "GSEA Measurement")

ggsave("hallmarks.png")

library(pheatmap)

means <- as.data.frame(rowMeans(hallmark.ssGSEA))
colnames(means) <- c("Mean ssGSEA Value")
means$id<-rownames(means)
hallmark.ssGSEA$id <- rownames(hallmark.ssGSEA)
foo <- full_join(hallmark.ssGSEA, means)
rownames(foo) <- foo$id
foo <- select(foo, -id)
pheatmap(foo, color = viridis(100), cluster_cols = FALSE, fontsize_row = 7, fontsize_col = 10, border_color = NA, cellwidth = 7, cellheight = 7)
synStore(File('avg_hallmark_ssGSEA.txt', parentId = 'syn7818711'), used = c('syn5556216','syn6156140'), executed = this.file, activityName = 'ssGSEA')

