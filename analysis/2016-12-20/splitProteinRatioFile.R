library(synapseClient)
library(dplyr)

NF_all_sets_redo<-synGet(id='syn7349351')
filepath <- NF_all_sets_redo@filePath
NFredo <- read.table(file = filepath, sep = '\t', header = TRUE)

source.file <- NF_all_sets_redo@fileHandle$fileName
this.file <- 'https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2016-12-20/splitProteinRatioFile.R'

c <- make.names(colnames(NFredo))
colnames(NFredo) <- c

for(i in colnames(NFredo[,6:47])) {
  protein <- select(NFredo, Protein.Group:Sequence.Name)
  source <- c("originalFileName" = rep(source.file, nrow(protein)))
  protein2 <- c(i = NFredo[,i])
  protein <- cbind(protein, protein2, source)
  protein <- write.table(protein, file = paste('dermalNF_updated_proteomics_',i,'_normalized.txt', sep = ""), sep = "\t", row.names = FALSE)
  #synStore(File(path = paste('dermalNF_updated_proteomics_',x,'_normalized.txt', sep = ""), parentId = 'syn7349340'), used = 'syn7349351', executed = this.file)
}
  