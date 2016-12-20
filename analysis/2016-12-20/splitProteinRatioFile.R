source('../../bin/dermalNFData.R')
library(synapseClient)
library(dplyr)

NF_all_sets_redo<-synGet(id='syn7349351')
filepath <- NF_all_sets_redo@filePath
NFredo <- read.table(file = filepath, sep = '\t')

source.file <- NF_all_sets_redo@fileHandle$fileName
this.file = 

for(x in colnames(NFredo[,6:47])) {
  protein <- select(NFredo, V1:V5, x, V48)
  source <- c("originalFileName", rep(source.file, nrow(protein)-1))
  protein <- cbind(protein, source)
  protein <- write.table(file = paste('dermalNF_updated_proteomics_',x,'_normalized.txt'), sep = "\t")
  synStore(File(path = paste('dermalNF_updated_proteomics_',x,'_normalized.txt'), parentID = 'syn7349340'), used)
}
