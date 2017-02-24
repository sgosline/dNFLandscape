source("../../bin/wgsAnalysis.R")
library(dplyr)
this.file = "https://raw.githubusercontent.com/allaway/dNFLandscape/776770a72dba97d0e7e20e97dc47d56d7b48b3df/analysis/2017-02-23/TestForMutantDrugEnrichment.R"

testgenes<-filter(all.gene.muts, Effect %in% c("MODERATE", "HIGH"))
testgenes<-count(dplyr::select(distinct(dplyr::select(testgenes, Sample, Gene)), Gene))
factor(testgenes$Gene, levels = testgenes$Gene[order(desc(testgenes$freq))])

list<-testgenes$Gene[order(desc(testgenes$freq))]

## pull drug data and filter for human targets, and eliminate drugs with
## 0 quantitative effects measured
drugdat <- synTableQuery("SELECT * FROM syn7341038")
drugdat <- as.data.frame(drugdat@values)
drugdat.filtered <- filter(drugdat, Organism == "Homo sapiens")
drugdat.filtered <- filter(drugdat.filtered, N_quantitative != 0)

## obtain data to map Uniprot id with Hugo Genes
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
bm <-
  getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"),
        mart = mart)
bm <- filter(bm, uniprot_swissprot != "")
colnames(bm) <- c("Hugo_Gene", "Uniprot_accession_numbers")

## map uniprot targets to hugo genes - generate list of dataframes, each
## data frame lists targets of drug
drugdat.filtered$Uniprot_accession_numbers <-
  sub(",.", "", drugdat.filtered$Uniprot_accession_numbers)
drugdat.filtered <-
  left_join(drugdat.filtered, bm, by = "Uniprot_accession_numbers")
listofdrugtargets <-
  dlply(.data = drugdat.filtered, .variables = "Structure_ID")

## some uniprot ids do not successfully map to hugo genes, note them
## here, these will not show up in analysis
unannotated <-
  filter(drugdat.filtered, is.na(drugdat.filtered$Hugo_Gene))

## pull evotec data to map compound names to structure IDs
compound.data <- synTableQuery("SELECT * FROM syn8118065")
compound.data <- as.data.frame(compound.data@values)
compound.data <-
  compound.data[!duplicated(compound.data$Structure_ID),]

## make more stringent version of above drug data (increased evidence
## for hitting target, more potency towards target) forget everything
## requiring more than 50uM drug to have an effect, would be very
## challenging to deliver that concentration of drug to target in tumor
drugdat.stringent <- filter(drugdat.filtered, N_quantitative > 1)
drugdat.stringent <-
  filter(drugdat.filtered, MinActivity_nM < 50000)
listofdrugtargets.stri <-
  dlply(.data = drugdat.stringent, .variables = "Structure_ID")

print(detectCores())

## core function to test for enriched drug targets (hypergeometric test)
TestForDrugTargets <- function(comut) {
  allcomuts <- unique(comut)
  hyper <- mclapply(listofdrugtargets.stri, function(x) {
    N <- length(allcomuts)
    B <- nrow(x$Hugo_Gene)
    lambdas <- as.integer((allcomuts %in% x$Hugo_Gene))
    mHG <- mHG.test(lambdas)$p.value
  }, mc.cores=detectCores())
  
  Structure_ID <- names(listofdrugtargets.stri)
  hypergeo_pval <- t(bind_rows(hyper))
  hyper.df <- as.data.frame(cbind(Structure_ID, hypergeo_pval))
  names(hyper.df) <- c("Structure_ID", "Hypergeo_pval")
  compound.data$Structure_ID <-
    as.character(compound.data$Structure_ID)
  hyper.annot <-
    left_join(hyper.df, compound.data, by = "Structure_ID")
}

hyper<-TestForDrugTargets(list)
hyper<-as.data.frame(hyper)

write.table(hyper,'cNF_somatic_muts_drug_enrichment.txt', sep = "\t")
synStore(File('cNF_somatic_muts_drug_enrichment.txt', parentId = "syn8319243"), used = c("syn7341038", "syn8118065"), executed = this.file)

library(synapseClient)
library(ggplot2)
library(rJava)
library(rcdk)
synapseLogin()
library(pheatmap)
library(viridis)

hyper<-read.table(synGet("syn8319288")@filePath, header = TRUE, sep = "\t")

###NF1 cNF somatic mutation hits 
sigs<-filter(hyper, Hypergeo_pval<0.05)
sig.data<-dplyr::select(sigs, Structure_ID, Original_molecule_SMILES, Supplier_Data_1,Supplier_Data_2,Supplier_Data_3)

mol<-parse.smiles(as.character(sig.data$Original_molecule_SMILES))

#write.molecules(mol, filename = 'NF1_TCGA_drugs.sdf')
#view.molecule.2d(mol)

fps <- lapply(mol, get.fingerprint, type="extended")
fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
rownames(fp.sim) <- sig.data$Structure_ID
colnames(fp.sim) <- sig.data$Structure_ID
fp.dist <- 1- fp.sim
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000))

clust<-hclust(as.dist(fp.dist))
plot(clust)


###NF1 TCGA Comutant Hits
x<-synGet("syn8292948")@filePath
data<-read.table(x)
sigs<-filter(data, Hypergeo_pval<=0.0005)
sig.data<-select(sigs, Structure_ID, Original_molecule_SMILES, Supplier_Data_1,Supplier_Data_2,Supplier_Data_3)
sig.data<-distinct(sig.data)

mol<-parse.smiles(as.character(sig.data$Original_molecule_SMILES))

#write.molecules(mol, filename = 'NF1_TCGA_drugs.sdf')
#view.molecule.2d(mol)

fps <- lapply(mol, get.fingerprint, type="extended")
fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
rownames(fp.sim) <- sig.data$Structure_ID
colnames(fp.sim) <- sig.data$Structure_ID
fp.dist <- 1- fp.sim

png("tanimotoSimilarityOfHits.png")
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000))
dev.off()

synStore(File('tanimotoSimilarityOfHits.png', parentId = "syn8319243"), used = c("syn7341038", "syn8118065", "syn8292948"), executed = this.file)



