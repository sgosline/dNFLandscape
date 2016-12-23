source("../../bin/wgsAnalysis.R")
source("../../analysis/2016-09-20/geneSampleMatrix.R")

som.vars<-somaticGeneSampleMatrix()
germ.vars<-germlineGeneSampleMatrix()

som.vars<-as.matrix(som.vars)

true.mutations<-list()

for(i in colnames(som.vars)){
 mutation<-which(som.vars[,i], arr.ind = TRUE)
 true.mutations[[i]] <- names(mutation)
 print(names(mutation))
 print(i)
}

truemut<-plyr::ldply(true.mutations, rbind)
write.table(truemut, file="truemut.m2", sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")

germ.vars <- as.matrix(germ.vars)

true.mutations<-list()

for(i in colnames(germ.vars)){
  mutation<-which(germ.vars[,i], arr.ind = TRUE)
  true.mutations[[i]] <- names(mutation)
  print(names(mutation))
  print(i)
}

truemut2<-plyr::ldply(true.mutations, rbind)
write.table(truemut2, file="truegermmut.m2", sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")

merge<-dplyr::full_join(truemut,truemut2, by = ".id")
write.table(merge, file="allmut.m2", sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE, na = "")
