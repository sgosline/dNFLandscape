##import dNF and ENCODE skin sample gene expression data 
source("../2016-11-16/dNF_ENCODE_RNAseq.R")

##all gene sets pulled from Gene Ontology

##first plot heatmap of genes related to mast cell chemotaxis
##import gene set and filter gene expression list to only include these data
mast.chemotaxis <- c("VEGFD", "VEGFC", "PIK3CD", "CCL11", "RAC1", "SWAP70", "PGF", "KIT", "CHGA") ##http://amigo.geneontology.org/amigo/term/GO:0002551
mast.chemotax.df <- dplyr::filter(norm.df, Gene %in% mast.chemotaxis)
mast.chemotax.df <- tidyr::spread(mast.chemotax.df, Gene, Expression)

##set gene names as rownames rather than first row, convert from dataframe to matrix
mast.chemotax.col <- mast.chemotax.df[,-1]
rownames(mast.chemotax.col) <- mast.chemotax.df[,1]
mast.chemotax.m <- as.matrix(mast.chemotax.col)

##make a dataframe to label all tumors with correct size measurement in mm, this is for annotation. set ENCODE samples to 0 as they are not tumors
samples <- rownames(mast.chemotax.m)
sizes <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,13,14,5,7,11,16,20,15,8,5,10,7,8,6,15,10,15,13,10,10,13,10,12,5,11,20,10,3,18,12,20,13)
tumorsize.m <- matrix(data = sizes, nrow = 66, ncol = 1, dimnames = (list(samples, "Size")))
tumorsize.df <- as.data.frame(tumorsize.m)

samples <- rownames(mast.chemotax.m)

##print pdf to file of heatmap
pdf("mast.chemotax.pdf")
pheatmap::pheatmap(mast.chemotax.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)
dev.off()

##repeat process for genes invovled in mast cell proliferation
mast.prolif <- c("LYN", "IL3", "KITLG", "KIT") ##http://amigo.geneontology.org/amigo/term/GO:0002551
mast.prolif.df <- dplyr::filter(norm.df, Gene %in% mast.prolif)
mast.prolif.df <- tidyr::spread(mast.prolif.df, Gene, Expression)

mast.prolif.col <- mast.prolif.df[,-1]
rownames(mast.prolif.col) <- mast.prolif.df[,1]
mast.prolif.m <- as.matrix(mast.prolif.col)

samples <- rownames(mast.prolif.m)

pdf("mast.prolif.pdf")
pheatmap::pheatmap(mast.prolif.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)
dev.off()

##repeat process for genes involved in mast cell mediated immune responses
mast.med.immunity <- c("SERPINB9", "GAB2", "IL4", "FCER1G", "RAC2", "FOXF1", "GATA2", "LYN", "PLA2G3", "IL13", "HMOX1", "C12orf4", "IL13RA2", "IL4R", "LAT", "LAT2", "SPON2", 
                       "ADORA2B", "MRGPRX2", "MILR1", "PDPK1", "FGR", "NR4A3", "ZAP70", "UNC13D", "FER", "PIK3CD", "LGALS9", "SNAP23", "STXBP2", "PIK3CG", "CPLX2", "SYK", "VAMP7",
                       "S100A13", "VAMP8", "CD300A", "FES", "FCER1A", "RABGEF1", "KIT", "BTK", "CD84", "RASGRP1", "SNX4") ##http://amigo.geneontology.org/amigo/term/GO:0002448
mast.med.immunity.df <- dplyr::filter(norm.df, Gene %in% mast.med.immunity)
mast.med.immunity.df <- tidyr::spread(mast.med.immunity.df, Gene, Expression)

mast.med.immunity.col <- mast.med.immunity.df[,-1]
rownames(mast.med.immunity.col) <- mast.med.immunity.df[,1]
mast.med.immunity.m <- as.matrix(mast.med.immunity.col)

samples <- rownames(mast.med.immunity.m)

pdf("mast.med.immunity.pdf")
pheatmap::pheatmap(mast.med.immunity.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)
dev.off()

##repeat process for genes involved in negative regulation of mast cell activation
neg.reg.mast.activation <- c("FOXF1", "HMOX1", "IL13RA3", "FER", "LGALS9", "CD300A", "RABGEF1", "CD84") ##http://amigo.geneontology.org/amigo/term/GO:0033004
neg.reg.mast.activation.df <- dplyr::filter(norm.df, Gene %in% neg.reg.mast.activation)
neg.reg.mast.activation.df <- tidyr::spread(neg.reg.mast.activation.df, Gene, Expression)

neg.reg.mast.activation.col <- neg.reg.mast.activation.df[,-1]
rownames(neg.reg.mast.activation.col) <- neg.reg.mast.activation.df[,1]
neg.reg.mast.activation.m <- as.matrix(neg.reg.mast.activation.col)

samples <- rownames(neg.reg.mast.activation.m)

pdf("neg.reg.mast.activation.pdf")
pheatmap::pheatmap(neg.reg.mast.activation.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)
dev.off()

##repeat process for genes involved in positive regulation of mast cell activation
pos.reg.mast.activation <- c("GAB2", "IL4", "FCER1G", "GATA2", "IL13", "IL4R", "ADORA2B", "FGR", "ZAP70", "SYK") ##http://amigo.geneontology.org/amigo/term/GO:0033005
pos.reg.mast.activation.df <- dplyr::filter(norm.df, Gene %in% pos.reg.mast.activation)
pos.reg.mast.activation.df <- tidyr::spread(pos.reg.mast.activation.df, Gene, Expression)

pos.reg.mast.activation.col <- pos.reg.mast.activation.df[,-1]
rownames(pos.reg.mast.activation.col) <- pos.reg.mast.activation.df[,1]
pos.reg.mast.activation.m <- as.matrix(pos.reg.mast.activation.col)

samples <- rownames(pos.reg.mast.activation.m)

pdf("pos.reg.mast.activation.pdf")
pheatmap::pheatmap(pos.reg.mast.activation.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)
dev.off()

##combine heatmaps for genes causing positive and negative activation of mast cells
mast.activation <- dplyr::bind_cols(pos.reg.mast.activation.df, neg.reg.mast.activation.df[2:8])

##create data frame to enable annotation of genes by direction of mast cell activation (+ = activator, - = repressor)
mast.activation.genes <- colnames(mast.activation[2:18])
mast.activation.direction <- c("+", "+", "+", "+", "+", "+", "+", "+", "+", "+", "-", "-", "-", "-", "-", "-", "-")
mast.activation.annotation <- data.frame(mast.activation.direction, row.names = mast.activation.genes)

mast.activation.col <- mast.activation[,-1]
rownames(mast.activation.col) <- mast.activation[,1]
mast.activation.m <- as.matrix(mast.activation.col)

##generate heatmap, preventing column clustering to so that + genes and - genes are segregated
samples <- rownames(mast.activation.m)

pdf("reg.mast.activation.pdf")
pheatmap::pheatmap(mast.activation.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", 
                   annotation_row = tumorsize.df, annotation_col = mast.activation.annotation, cluster_cols = FALSE)
dev.off()

##repeat previous heatmap with genes regulating mast cell apoptosis
pos.reg.mast.apop <- c("FCER1G", "ADAM17", "KITLG") ##http://amigo.geneontology.org/amigo/term/GO:0033024
pos.reg.mast.apop.df <- dplyr::filter(norm.df, Gene %in% pos.reg.mast.apop)
pos.reg.mast.apop.df <- tidyr::spread(pos.reg.mast.apop.df, Gene, Expression)

pos.reg.mast.apop.col <- pos.reg.mast.apop.df[,-1]
rownames(pos.reg.mast.apop.col) <- pos.reg.mast.apop.df[,1]
pos.reg.mast.apop.m <- as.matrix(pos.reg.mast.apop.col)

samples <- rownames(pos.reg.mast.apop.m)

pdf("pos.reg.mast.apop.pdf")
pheatmap::pheatmap(pos.reg.mast.apop.m, cellwidth = 5, cellheight = 5, fontsize = 5, treeheight_col= 40, treeheight_row = 40, border_color = "grey", annotation_row = tumorsize.df)
dev.off()
