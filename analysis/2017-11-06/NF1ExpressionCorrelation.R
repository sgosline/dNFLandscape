source("../../dermalNF/bin/dermalNFData.R")
library(plyr)

this.file = "https://github.com/allaway/dNFLandscape/blob/master/analysis/2017-11-06/NF1ExpressionCorrelation.R"

res = rna_count_matrix(stored=TRUE,doNorm=FALSE,minCount=0,doLogNorm=FALSE,doVoomNorm=FALSE)

moreids<-synTableQuery("SELECT * FROM syn5556216")@values %>% filter(!is.na(RNASeq)) %>% 
  dplyr::select(sampleIdentifier, RNASeq) %>% 
  arrange(RNASeq)

colnames(res) <- gsub("X", "", colnames(res))
colnames(res) <- as.character(gsub("\\.", "-", colnames(res)))
colnames(res) <- moreids$sampleIdentifier[colnames(res) == moreids$RNASeq]

moreids$patient <- gsub("tumor\\d+$", "", moreids$sampleIdentifier)

corres <- as.data.frame(t(res))

cors <- lapply(colnames(corres), function(x) {
  foo <- cor.test(corres[[x]], corres$NF1)
  c(foo$estimate, foo$p.value)
})

names(cors) <- colnames(corres)
cors <- ldply(cors)
cors$bh <- p.adjust(cors$V1, method = "BH")
write.table(cors, "cNF-NF1-Correlation.txt", sep = "\t")
synStore(File("cNF-NF1-Correlation.txt", parentId = "syn11664760"), 
         used = c("syn5556216","syn5051784"), executed = this.file)

##top 200 genes
up <- cors %>% filter(bh <0.05) %>% top_n(200, cor)
down <- cors %>% filter(bh < 0.05) %>% top_n(200, -cor)

library(enrichR)
dbs <- c("GO_Biological_Process_2017b",
         "GO_Molecular_Function_2017b",
         "GO_Cellular_Component_2017b",
         "KEGG_2016",
         "Human_Phenotype_Ontology")

cnf.up <- enrichr(up$.id, databases = dbs)
cnf.down <- enrichr(down$.id, databases = dbs)


cnf.up <- ldply(cnf.up)

cnf.up$SigTerm <- " "
cnf.up$SigTerm[cnf.up$Adjusted.P.value<0.05] <- cnf.up$Term[cnf.up$Adjusted.P.value<0.05]

library(ggplot2)
library(ggrepel)

ggplot(cnf.up %>% filter(.id != "KEGG_2016"), aes(x = Z.score, y = -log(Adjusted.P.value))) + 
  geom_point(aes(color = .id)) +
  geom_text_repel(data = cnf.up %>% filter(Adjusted.P.value < 0.05, .id != "KEGG_2016"), aes(label = SigTerm))

cnf.down <- ldply(cnf.down)

cnf.down$SigTerm <- " "
cnf.down$SigTerm[cnf.down$Adjusted.P.value<0.05] <- cnf.down$Term[cnf.down$Adjusted.P.value<0.05]

ggplot(cnf.down %>% filter(.id != "KEGG_2016"), aes(x = Z.score, y = -log(Adjusted.P.value))) + 
  geom_point(aes(color = .id)) +
  geom_text_repel(data = cnf.down %>% filter(Adjusted.P.value < 0.05, .id != "KEGG_2016"), 
                  aes(label = SigTerm),
                  size = 2.5)


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


for(i in dbs){

  cnf.up <- enrichr(up$.id, databases = i)
  cnf.down <- enrichr(down$.id, databases = i)
  
  cnf.up <- ldply(cnf.up)
  cnf.down <- ldply(cnf.down)
  
  if(nrow(cnf.up)>1){
    cnf.up$SigTerm <- " "
    cnf.up$SigTerm[cnf.up$Adjusted.P.value<0.05] <- cnf.up$Term[cnf.up$Adjusted.P.value<0.05]
    
    ggplot(cnf.up , aes(x = Z.score, y = -log(Adjusted.P.value))) +
    geom_point(aes(color = .id, size = Combined.Score/2)) +
      geom_text_repel(data = cnf.up %>% filter(Adjusted.P.value < 0.05), aes(label = SigTerm))
    
    ggsave(paste0(i,"_enrichr_positivelycorrelatedgenes.png"), height = 8, width = 10)
    synStore(File(paste0(i,"_enrichr_positivelycorrelatedgenes.png"), parentId = "syn11664964"),
             used = c("syn5556216","syn11664965"), executed = this.file)
  }
  
  if(nrow(cnf.down)>1){  
    cnf.down$SigTerm <- " "
    cnf.down$SigTerm[cnf.down$Adjusted.P.value<0.05] <- cnf.down$Term[cnf.down$Adjusted.P.value<0.05]
    
    ggplot(cnf.down, aes(x = Z.score, y = -log(Adjusted.P.value))) + 
      geom_point(aes(color = .id, size = Combined.Score/2)) +
      geom_text_repel(data = cnf.down %>% filter(Adjusted.P.value < 0.05), 
                      aes(label = SigTerm),
                      size = 2.5)
    
    ggsave(paste0(i,"_enrichr_negativelycorrelatedgenes.png"), height = 8, width = 10)
    synStore(File(paste0(i,"_enrichr_negativelycorrelatedgenes.png"), parentId = "syn11664964"),
             used = c("syn5556216","syn11664965"), executed = this.file)
  }
}


###PCSF correlated genes 

library(PCSF)

data("STRING")
ppi <- construct_interactome(STRING)

sigcors <- cors %>% filter(bh < 0.05, cor < -0.5)

terminals <- rank(-sigcors$cor)
names(terminals) <- sigcors$.id

subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.05)
plot(subnet)

subnet <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 2, b = 1, mu = 0.0005) ##bit of hairball
subnet <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 2, b = 1, mu = 0.005)

plot(subnet) %>% visSave("PCSF_NF1_negatively_correlated_genes.html", selfcontained = T)

synStore(File("PCSF_NF1_negatively_correlated_genes.html", parentId = "syn11664974"),
         used = c("syn11664965"), executed = this.file)

res <- enrichment_analysis(subnet)
plot(res$subnet) %>% visSave("PCSF_enrichr_Annotated.html", selfcontained = T)

synStore(File("PCSF_enrichr_Annotated.html", parentId = "syn11664974"),
         used = c("syn11664965"), executed = this.file)
