library(plyr)
library(dplyr)
library(tidyr)
library(tibble)
library(PCSF)
library(synapseClient)
synapseLogin()

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/NF_Cancer_Landscape/master/analysis/#######"

oncotree <- read.table(synGet("syn7392892", version = 117)@filePath, sep = "\t", header = T, quote = "\"", comment.char = "") %>% 
  select("SAMPLE_ID","ONCOTREE_CODE")

genie <- read.table(synGet("syn11316113")@filePath, header = T, row.names = 1) %>% t() %>% as.data.frame() %>% rownames_to_column("gene")
genie.nf1 <- genie %>% filter(gene == "NF1")  %>% t() %>% as.data.frame() %>% rownames_to_column("patient")
genie.nf1 <- genie.nf1[-1,] %>% filter(V1 == 1)

genie <- genie %>% 
  select(gene, one_of(genie.nf1$patient)) %>% 
  column_to_rownames("gene") %>% 
  rownames_to_column("gene") %>% 
  gather("SAMPLE_ID", "mutation", -gene) %>% 
  left_join(oncotree)

genie.sum <- genie %>% 
  group_by(gene, ONCOTREE_CODE) %>% 
  dplyr::summarise(mutation2 = sum(mutation, na.rm = T)) %>% 
  filter(ONCOTREE_CODE != "NA") %>% 
  ungroup()

data("STRING")
ppi <- construct_interactome(STRING)

####################################### GBM NF1 Mutant Samples
genie.gbm <- genie.sum %>% filter(ONCOTREE_CODE == "GBM", mutation2 > 0)

terminals <- genie.gbm$mutation2
names(terminals) <- genie.gbm$gene

# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 2, b = 1, mu = 0.0005) ##bit of hairball
# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 5, b = 1, mu = 0.0005) ##slighly less of hairball
# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 5, b = 10, mu = 0.0005) ##no obvious difference
# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 5, b = 1000, mu = 0.0005) ##slighty longer branches but no huge change
# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1000, b = 1, mu = 0.0005) ##more branches
# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1000, b = 1, mu = 0.005) ##lose most of steiner nodes
# subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1000, b = 1, mu = 0.0000005) ##UBC appears
#subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 0.1, b = 1, mu = 0.0005) ## no steiner nodes
subnet.gbm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1, b = 1, mu = 0.0005) ## not super hairbally

plot(subnet.gbm)

res.gbm <- enrichment_analysis(subnet.gbm)
plot(res.gbm$subnet)

####################################### LGG NF1 Mutant Samples
genie.lgg <- genie.sum %>% filter(ONCOTREE_CODE == "LGGNOS", mutation2 > 0)

terminals <- genie.lgg$mutation2
names(terminals) <- genie.lgg$gene

subnet.lgg <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1, b = 1, mu = 0.0005) 

plot(subnet.lgg)
res.lgg <- enrichment_analysis(subnet.lgg)
plot(res.lgg$subnet)

####################################### SKCM NF1 Mutant Samples
genie.skcm <- genie.sum %>% filter(ONCOTREE_CODE == "SKCM", mutation2 > 0)

terminals <- genie.skcm$mutation2
names(terminals) <- genie.skcm$gene

subnet.skcm <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1, b = 1, mu = 0.0005) ##bit of hairball

plot(subnet.skcm)
res.skcm <- enrichment_analysis(subnet.skcm)
plot(res.skcm$subnet)

####################################### LUAD NF1 Mutant Samples
genie.luad <- genie.sum %>% filter(ONCOTREE_CODE == "LUAD", mutation2 > 0)

terminals <- genie.luad$mutation2
names(terminals) <- genie.luad$gene

subnet.luad <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1, b = 1, mu = 0.0005) ##bit of hairball


plot(subnet.luad)
res.luad <- enrichment_analysis(subnet.luad)
plot(res.luad$subnet)


####################################### UEC NF1 Mutant Samples
genie.uec <- genie.sum %>% filter(ONCOTREE_CODE == "UEC", mutation2 > 0)

terminals <- genie.uec$mutation2
names(terminals) <- genie.uec$gene

subnet.uec <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 1, b = 1, mu = 0.0005) ##bit of hairball


plot(subnet.uec)
res.uec <- enrichment_analysis(subnet.uec)
plot(res.uec$subnet)

####################################### PanCan

genie <- genie.sum %>% group_by(gene) %>% dplyr::summarise(mutationall = sum(mutation2, na.rm = T)) %>% ungroup() %>% filter(mutationall > 0)

terminals <- genie$mutationall
names(terminals) <- genie$gene

subnet <- PCSF_rand(ppi, terminals, n = 10, r = 1, w = 10, b = 1, mu = 0.005)

plot(subnet)
plot(subnet2)

res <- enrichment_analysis(subnet)
plot(res$subnet)

## take all genes that are mutated in non nf1 mutant cancer
## 100 networks 
## weight each node (number of times mutated/all mutants??)
##hotnet statistics to determine if network meaningful 
## more forests less trees 
