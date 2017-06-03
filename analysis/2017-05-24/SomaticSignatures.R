 library(dplyr)
library(tidyr)
library(stringr)
library(synapseClient)
library(SomaticSignatures)
synapseLogin()

vcf <- readVcfAsVRanges(synGet("syn5555584")@filePath)

sca_vr <- VRanges(
  ,
  ranges = ir,
  ref = sca_data$Reference_Allele,
  alt = sca_data$Variant_Allele,
  sampleNames = sca_data$Tumor_Sample_Barcode),
  seqinfo = rep("cNFgermline", 20039),
  study = rep("cNFgermline", 20039))

VRanges(
  seqnames = chr
)
