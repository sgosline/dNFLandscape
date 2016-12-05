##reformatting data 
mast.up <- unique(read.table("Ito_2012_MAST_UP.txt", header = TRUE))
mast.up <- toupper(mast.up[,1])
mast.up <- as.data.frame(mast.up)
colnames(mast.up) <- c("Gene")
mast.down <- unique(read.table("Ito_2012_MAST_DOWN.txt", header = TRUE))
mast.down <- toupper(mast.down[,1])
mast.down <- as.data.frame(mast.down)
colnames(mast.down) <- c("Gene")
direction <- c(rep.int("UP", 151), rep.int("DOWN", 346))
direction <- as.data.frame(direction)
mast.activation <- dplyr::bind_rows(mast.up, mast.down)
mast.activation <- dplyr::bind_cols(mast.activation, direction)


##get gene expression data 
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

##find overlap in mast cell activation set and gene expression data
mast.activation.Itosig <- dplyr::semi_join(comb.norm, mast.activation, by = "Gene")
mast.activation.labels <- dplyr::semi_join (mast.activation, mast.activation.Itosig, by = "Gene")

rownames(mast.activation.labels) <- mast.activation.labels$Gene
dirnum <- c(rep.int(1, 116), rep.int(0, 256))
mast.activation.labels <- cbind(mast.activation.labels, dirnum)
mast.activation.labels <- dplyr::select(mast.activation.labels, dirnum)


rownames(mast.activation.Itosig) <- mast.activation.Itosig$Gene
mast.activation.Itosig <- mast.activation.Itosig[,1:66]

pdf("mast.activation.Itosig.pdf", height = 15)
pheatmap::pheatmap(mast.activation.Itosig, cluster_rows = FALSE, cellwidth = 5, cellheight = 2, fontsize = 2, 
   treeheight_col= 40, treeheight_row = 40, border_color = "grey" , annotation_row = mast.activation.labels)
pheatmap::pheatmap(mast.activation.Itosig, cellwidth = 5, cellheight = 2, fontsize = 2, 
                   treeheight_col= 40, treeheight_row = 40, border_color = "grey" , annotation_row = mast.activation.labels)
pheatmap::pheatmap(mast.activation.Itosig[1:116,1:66], cellwidth = 5, cellheight = 2, fontsize = 2, 
                   treeheight_col= 40, treeheight_row = 40, border_color = "grey")
pheatmap::pheatmap(mast.activation.Itosig[117:372,1:66], cellwidth = 5, cellheight = 2, fontsize = 2, 
                   treeheight_col= 40, treeheight_row = 40, border_color = "grey")
pheatmap::pheatmap(mast.activation.Itosig[,1:33], cellwidth = 5, cellheight = 2, fontsize = 2, cluster_rows = FALSE, 
                   treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = mast.activation.labels)
pheatmap::pheatmap(mast.activation.Itosig[,1:33], cellwidth = 5, cellheight = 2, fontsize = 2, cluster_rows = TRUE, 
                   treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = mast.activation.labels)

dev.off()

