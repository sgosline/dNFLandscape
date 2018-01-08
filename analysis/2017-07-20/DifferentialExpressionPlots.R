library(synapseClient)
synapseLogin()
library(tidyverse)
library(ggplot2)
library(enrichR)

degs <- read.table(synGet("syn11395910")@filePath, sep = "\t", header = T)
degs.sig <- degs %>% filter(padj < 0.1)

degs.05 <- degs %>% filter(padj < 0.05)
degs.sum <- as.data.frame(table(degs.05$gene)) %>% arrange(Freq)
degs.sum$Var1 <- factor(degs.sum$Var1, levels = degs.sum$Var1[order(-degs.sum$Freq)])

ggplot(degs.sum %>% filter(Freq > 1)) +
  geom_bar(aes(x = Var1, y = Freq, fill = Freq), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

enrich.res <- as.data.frame(c())
dbs <- c("GO_Cellular_Component_2017b",
         "GO_Molecular_Function_2017b",
         "GO_Biological_Process_2017b",
         "PPI_Hub_Proteins",
         "KEGG_2016",
         "CHEA_2013",
         "CORUM",
         "Disease_Signatures_from_GEO_down_2014",
         "Disease_Signatures_from_GEO_up_2014",
         "HumanCyc_2015",
         "MGI_Mammalian_Phenotype_Level_3",
         "MGI_Mammalian_Phenotype_Level_4",
         "MSigDB_Computational",
         "MSigDB_Oncogenic_Signatures",
         "OMIM_Disease",
         "Reactome_2013",
         "WikiPathways_2013")

for(i in unique(degs.sig$test)){
  foo <- degs.sig %>% filter(test == i)
  res <- enrichr(foo$gene, databases = dbs)
  res <- ldply(res)
  res$test <- i
  enrich.res <- bind_rows(enrich.res, res)
}

enrich.tab <- enrich.res %>% filter(Adjusted.P.value < 0.1) ### nothing overlaps at p <0.05
enrich.tab <- as.data.frame(table(enrich.tab$Term)) %>% arrange(Freq)
enrich.tab$Var1 <- factor(enrich.tab$Var1, levels = enrich.tab$Var1[order(-enrich.tab$Freq)])

ggplot(enrich.tab %>% filter(Freq > 1)) +
  geom_bar(aes(x = Var1, y = Freq, fill = Freq), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
