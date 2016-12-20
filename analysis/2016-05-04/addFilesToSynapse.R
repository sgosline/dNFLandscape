library(synapseClient)
synapseLogin()

diffexScriptlist=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-05-04/encodeSkinAnalysis_ExprOnly.R'),
  list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/bin/encodeSkinRNASeq.R'),list(url='https://raw.githubusercontent.com/Sage-Bionetworks/dermalNF/master/bin/dermalNFData.R'))

txtfiles=list.files('.')
txtfiles=txtfiles[grep("txt",txtfiles)]

for(file in txtfiles){
    synStore(File(file,parentId='syn7518454'),executed=diffexScriptlist,used=list(list(entity='syn6035999'),list(entity='syn5579598')))
}
