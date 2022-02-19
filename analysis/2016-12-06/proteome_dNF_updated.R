library(synapseClient)
source('../../dermalNF/bin/dermalNFData.R')

#get updated proteomics data
syn7349351 <- synGet(id='syn7349351')

annotations <- protein_annotations()

protein.unnormalized <- prot_unnormalized()

