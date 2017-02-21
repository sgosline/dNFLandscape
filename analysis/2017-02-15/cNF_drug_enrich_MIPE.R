library(synapseClient)
library(plyr)
library(dplyr)
library(mHG)
library(ggplot2)

## pull genie NF1 comutation data from synapse, break into list of dfs
## by cancer type
topNF1comut.genie <- as.data.frame(synTableQuery("SELECT * FROM syn8073830")@values)
listofcomuts <- dlply(.data = topNF1comut.genie, .variables = "ONCOTREE_CODE")

topNF1comut.TCGA <- as.data.frame(synTableQuery("SELECT * FROM syn8069190")@values)
listofcomuts.TCGA <- dlply(.data = topNF1comut.TCGA, .variables = " Cancer_Type")

cNF_genes <- c("CYLD", "FLI1", "HOXA13", "MN1", "NSD1", "RHOA", "EXT1", "IL6ST", "MAF",
               "NFKB2", "TSC1", "COL1A1", "DDX5", "RAD21", "ELL", "PBRM1", "RUNX1", "NF1",
               "MNX1")
mut_freq <- c(2,2,2,2,2,2,3,3,3,3,3,5,5,7,9,9,12,17,22)
cNF_genes <- as.data.frame(cNF_genes)
mut_freq <- as.data.frame(mut_freq)
cNF<-cbind(cNF_genes, mut_freq)

## pull evotec data to map compound names to structure IDs
druglist <- read.table(file=synGet("syn5522649")@filePath, sep = ",", header = TRUE, fill = TRUE, quote = "\"")
compound.data <- select(druglist, name, target)
compound.data <- filter(compound.data, !is.na(target))
listofdrugtargets.stri <- dlply(compound.data, .variables = "name")

## core function to test for enriched drug targets (hypergeometric test)
TestForDrugTargets <- function(comut) {
  allcomuts <- dplyr::select(comut, Hugo_Symbol, mut_freq_in_cancer_type)
  allcomuts <- arrange(allcomuts, desc(mut_freq_in_cancer_type))
  allcomuts <- unique(allcomuts)
  hyper <- lapply(listofdrugtargets.stri, function(x) {
    N <- length(allcomuts)
    B <- nrow(x$Hugo_Gene)
    lambdas <- as.integer((allcomuts$Hugo_Symbol %in% x$Hugo_Gene))
    mHG <- mHG.test(lambdas)$p.value
  })
  
  name <- names(listofdrugtargets.stri)
  hypergeo_pval <- t(bind_rows(hyper))
  hyper.df <- as.data.frame(cbind(name, hypergeo_pval))
  names(hyper.df) <- c("name", "Hypergeo_pval")
  compound.data$name <- as.character(compound.data$name)
  hyper.annot <- left_join(hyper.df, compound.data, by = "name")
}

## lapply across all tumor types
hyper.MIPE.GENIE <- lapply(names(listofcomuts), function(x) {
  print(x)
  hyper.annot <- TestForDrugTargets(listofcomuts[[x]])
  cancer_type <- rep(x, nrow(hyper.annot))
  hyper.annot <- cbind(hyper.annot, cancer_type)
})

## consolidate back into one df, adjust pval for multi corrections, make
## df with significantly enriched compounds
hyper.MIPE.GENIE.df <- bind_rows(hyper.MIPE.GENIE)
hyper.MIPE.GENIE.df$pval_BHadj <- p.adjust(hyper.MIPE.GENIE.df$Hypergeo_pval, method = "BH")
sigs.mipe.genie <- filter(hyper.MIPE.GENIE.df, pval_BHadj < 0.01)
ids.g <- count(sigs.g$Structure_ID)

## now same for TCGA

## lapply across all tumor types
hyper.MIPE.TCGA <- lapply(names(listofcomuts.TCGA), function(x) {
  print(x)
  hyper.annot <- TestForDrugTargets(listofcomuts.TCGA[[x]])
  cancer_type <- rep(x, nrow(hyper.annot))
  hyper.annot <- cbind(hyper.annot, cancer_type)
})

## consolidate back into one df, adjust pval for multi corrections, make
## df with significantly enriched compounds
hyper.MIPE.TCGA.df <- bind_rows(hyper.MIPE.TCGA)
hyper.MIPE.TCGA.df$pval_BHadj <- p.adjust(hyper.MIPE.TCGA.df$Hypergeo_pval, method = "BH")
sigs.mipe.tcga <- filter(hyper.MIPE.TCGA.df, pval_BHadj < 0.05)
#write.table(sigs.tcga, file="NF1_coMutant_Targets_Enriched_Drugs.txt", sep = "\t")
#project <- synGet('syn7154892')
#tcresult <- as.tableColumns(sigs.tcga)
#cols <- tcresult$tableColumns
#schema <- TableSchema(name='NF1 co-mutant Targets - Enriched Drugs', columns=cols, parent=project)
#table <- Table(schema, as.integer(sigs.tcga)
#table <- synStore(table)

ids.tcga <- count(sigs.tcga, Structure_ID)
ids.tcga <- arrange(ids.tcga, desc(n))
ids.tcga$Structure_ID <- factor(ids.tcga$Structure_ID, levels = ids.tcga$Structure_ID[order(desc(ids.tcga$n))])
ggplot(data = ids.tcga[1:20,], aes(x = ids.tcga$Structure_ID[1:20], y = ids.tcga$n[1:20])) +
  geom_bar(stat = "identity")

colnames(cNF)<-c("Hugo_Symbol", "mut_freq_in_cancer_type")
hyper.cNF <- TestForDrugTargets(cNF)
cancer_type <- rep(x, nrow(hyper.annot))
hyper.annot <- cbind(hyper.annot, cancer_type)
