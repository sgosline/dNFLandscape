##goal is to format purity estimates in easy wa

library(synapseClient)
synapseLogin()
library(data.table)

estimate<-as.data.frame(t(read.table(synGet('syn5908274')@filePath,sep='\t')))
require(ggplot2)


patient<-synTableQuery("select Patient,Length_in_mm,RNASeq,TumorLocation from syn5556216")@values
estimate$Patient=as.factor(patient$Patient[match(rownames(estimate),patient$RNASeq)])
estimate$TumorSize=patient$Length_in_mm[match(rownames(estimate),patient$RNASeq)]
estimate$Location=patient$TumorLocation[match(rownames(estimate),patient$RNASeq)]

#ggplot(estimate)+geom_point(aes(x=StromalScore,y=ImmuneScore,col=TumorPurity))
#ggsave('ESTIMATEscoresAcrossSamples.png')

r=cor(estimate$ImmuneScore,estimate$StromalScore)

ggplot(estimate)+geom_point(aes(x=StromalScore,y=ImmuneScore,col=Patient,size=TumorPurity))+ggtitle(paste("Stromal-Immune R =",format(r,digits=3)))
ggsave('ESTIMATEscoresAcrossSamples.png')

r=cor(estimate$TumorPurity,estimate$TumorSize)
ggplot(estimate)+geom_point(aes(x=TumorPurity,y=TumorSize,col=Location))+ggtitle(paste("Purity-Size R =",format(r,digits=3)))
ggsave('ESTIMATEpurityscoresAcrossSamplesBySize.png')

write.table(estimate,file='estimateScoresWithOtherTumorData.txt',sep='\t',row.names=F,quote=F)
synStore(File('estimateScoresWithOtherTumorData.txt',parentId='syn5908270'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-05-05/formatTumorPurityEstimates.R')),used=list(list(entity='syn5908274')))

write.table(data.frame(Sample=rownames(estimate),Purity=estimate$TumorPurity),file='tumorPurityOnly.txt',row.names=F,quote=F,sep='\t')
synStore(File('tumorPurityOnly.txt',parentId='syn5908270'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-05-05/formatTumorPurityEstimates.R')),used=list(list(entity='syn5908274')))

synStore(File('ESTIMATEpurityscoresAcrossSamplesBySize.png',parentId='syn5908270'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-05-05/formatTumorPurityEstimates.R')),used=list(list(entity='syn5908274')))
synStore(File('ESTIMATEscoresAcrossSamples.png',parentId='syn5908270'),executed=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-05-05/formatTumorPurityEstimates.R')),used=list(list(entity='syn5908274')))