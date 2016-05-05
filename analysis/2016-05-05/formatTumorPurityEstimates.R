##goal is to format purity estimates in easy wa

library(synapseClient)
synapseLogin()
library(data.table)

estimate<-as.data.frame(t(read.table(synGet('syn5908274')@filePath,sep='\t')))
require(ggplot2)


patient<-synTableQuery("select Patient,Length_in_mm,RNASeq from syn5556216")@values
estimate$Patient=as.factor(patient$Patient[match(rownames(estimate),patient$RNASeq)])
estimate$TumorSize=patient$Length_in_mm[match(rownames(estimate),patient$RNASeq)]
#ggplot(estimate)+geom_point(aes(x=StromalScore,y=ImmuneScore,col=TumorPurity))
#ggsave('ESTIMATEscoresAcrossSamples.png')

ggplot(estimate)+geom_point(aes(x=StromalScore,y=ImmuneScore,col=Patient,size=TumorPurity))
ggsave('ESTIMATEscoresAcrossSamples.png')

r=cor(estimate$TumorPurity,estimate$TumorSize)
ggplot(estimate)+geom_point(aes(x=TumorPurity,y=TumorSize,col=Patient))+ggtitle(paste("R =",format(r,ndigits=4)))
ggsave('ESTIMATEscoresAcrossSamplesBySize.png')