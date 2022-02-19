library(synapseClient)
library(PCSF)
library(dplyr)
library(tidyr)
library(stringr)

synapseLogin()

this.file = "https://github.com/allaway/dNFLandscape/blob/master/analysis/###########"

mut<-read.table(synGet("syn5839666")@filePath, sep = "\t", header = TRUE, colClasses = "character", stringsAsFactors = TRUE)

mut$chr <- str_extract(mut$Chromosome, "\\d")

soms <- filter(mut, Mutation_Status == "Somatic", Mutation_Type != "Silent", Mutation_Type != "Intron") %>% 
  select(Hugo_Symbol, Sample_ID) %>% 
  distinct() %>% 
  group_by(Hugo_Symbol) %>% 
  dplyr::summarize(n = n()) %>% 
  ungroup() %>% 
  filter(n > 5)


data("STRING")
ppi <- construct_interactome(STRING)

terminals <- soms$n
names(terminals) <- soms$Hugo_Symbol

subnet <- PCSF(ppi, terminals, w = 2, b = 1, mu = 0.0005)
plot(subnet)

subnet <- PCSF_rand(ppi, terminals, n = 10, r = 0.1, w = 10, b = 5, mu = 0.05) #a bit hairball-y

plot(subnet)

res <- enrichment_analysis(subnet)
plot(res$subnet, node_size = 40)

