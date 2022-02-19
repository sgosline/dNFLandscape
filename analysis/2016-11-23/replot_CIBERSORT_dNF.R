##copy of 2016-09-16/rePlotCiberSort.R, adjusted for cibersort with dermals + skin

##replot cibersort
require(synapseClient)
synapseLogin()
library(dplyr)

cs.res<-read.delim(synGet('syn5809355')@filePath,sep=",",header=T, quote = "\"")
cs.res$Input.Sample <- c("3096-PBK-0001", "3096-PBK-0002", "3096-PBK-0004", "3096-PBK-0005", "3096-PBK-0006", "3096-PBK-0009",
                            "3096-PBK-0010", "3096-PBK-0011", "3096-PBK-0012", "3096-PBK-0013", "3096-PBK-0014", "3096-PBK-0015",
                            "3096-PBK-0016", "3096-PBK-0019", "3096-PBK-0020", "3096-PBK-0021", "3096-PBK-0022", "3096-PBK-0023", 
                            "3096-PBK-0024", "3096-PBK-0025", "3096-PBK-0026", "3096-PBK-0027", "3096-PBK-0028", "3096-PBK-0029", 
                            "3096-PBK-0032", "3096-PBK-0033", "3096-PBK-0034", "3096-PBK-0035", "3096-PBK-0039", "3096-PBK-0041",
                            "3096-PBK-0043", "3096-PBK-0122", "3096-PBK-0123")

##the PBK ids are missing from table, so need to query annotations
map<-read.table("3096-PBK-RNA Samples with RIN numbers.csv",sep=",",header=T)
res<-synQuery("select sampleIdentifier from entity where parentId=='syn5492805'")
names(res) <- c("sampleIdentifier", "SynapseID")
full.map<-left_join(res, map) %>% filter(Sample.Id != "NA")

full.map$Sample.Id <- gsub("-", ".", full.map$Sample.Id, fixed = TRUE)
full.map$Sample.Id <- paste("X",full.map$Sample.Id, sep = "")

full.map$sampleIdentifier <- gsub("patient", "Patient ", full.map$sampleIdentifier)
full.map$sampleIdentifier <- gsub("tumor", " Tumor ", full.map$sampleIdentifier)

full.map<-select(full.map, sampleIdentifier, Sample.Id)
colnames(cs.res)[1] <- "Sample.Id"
cs.res <- left_join(cs.res, full.map)
rownames(cs.res) <- cs.res$sampleIdentifier

this.file='https://github.com/allaway/dNFLandscape/blob/master/analysis/2016-11-23/replot_CIBERSORT_skin.R'

library(pheatmap)
pheatmap(cs.res[,2:23],annotation_row=cs.res[,24:26],cellheight=10,cellwidth = 10, filename = 'ciberSortRePlotted_dNF.png')
synStore(File('ciberSortRePlotted_dNF.png',parentId='syn5809348'), used =c('syn5492805', "syn5809355"), executed=list(list(url=this.file)))
