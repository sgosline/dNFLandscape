#Test pathway enrichment of various sets of mutations
source("../../bin/wgsAnalysis.R")

#2#evaluate most prevalent mutations at high impact, somatic/germline

##let's delineate all the tests we need to do
effects=c("HIGH","HIGH_MODERATE_LOW","HIGH_MODERATE")
write(as.character(unique(expr.gene.muts1$Gene)),file='allGeneBG.txt')

for(e in effects){
  p05_mutcounts<-getMutsAcrossGenes(expr.gene.muts05,effect=unlist(strsplit(e,split='_')),germLine=c("Germline"),
                                    som=c("LikelySomatic",'StrongSomatic'))
  
  p1_mutcounts<-getMutsAcrossGenes(expr.gene.muts1,effect=unlist(strsplit(e,split='_')),germLine=c("Germline"),
                                   som=c("LikelySomatic",'StrongSomatic'))
  
  write(as.character(p05_mutcounts$all$Gene),file=paste("allGenesWithp05Muts",e,'.txt',sep=''))
  write(as.character(p1_mutcounts$all$Gene),file=paste("allGenesWithp1Muts",e,'.txt',sep=''))
  
  ##which germline variants are shared across all samples? 
  com.gl.05<- as.character(p05_mutcounts$germline$Gene[which(p05_mutcounts$germline$nPatients==9)])
  
  com.gl.1<- as.character(p1_mutcounts$germline$Gene[which(p1_mutcounts$germline$nPatients==9)])
  
  ##which somatic variants are common across at least 2? 
  com.so.05<- as.character(p05_mutcounts$somatic$Gene[which(p05_mutcounts$somatic$nSamps>1)])
  com.so.1<- as.character(p1_mutcounts$somatic$Gene[which(p05_mutcounts$somatic$nSamps>1)])
  
  ##now wriet out all 4 gene lists
  f1=paste('genesWithGLVarsP05AcrossAllPats',e,'.txt',sep='')
  f2=paste('genesWithGLVarsP1AcrossAllPats',e,'.txt',sep='')
  write(com.gl.05,f1)
  write(com.gl.1,f2)
  
  f3=paste('genesWithSOVarsP05AcrossAtleast2Samps',e,'.txt',sep='')
  f4=paste('genesWithSOVarsP1AcrossAtleast2Samps',e,'.txt',sep='')
  write(com.so.05,f3)
  write(com.so.1,f4)
  
}
