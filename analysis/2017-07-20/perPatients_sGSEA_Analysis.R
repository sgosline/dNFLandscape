## this currently spits out GSVA not ssGSEA results 
source("../../bin/ssGSEAandGSVA_perPatient.R")
ids<-read.table(synGet("syn6156140")@filePath, sep = ",", header = T)

colnames(ids)[4] <- "RNASeq (Cufflinks)"
moreids<-synTableQuery("SELECT * FROM syn5556216")@values %>% filter(!is.na(RNASeq)) %>% 
  left_join(ids) %>% 
  dplyr::select(sampleIdentifier, Sample.Id) %>% 
  arrange(Sample.Id)

moreids$Sample.Id <- gsub("3096", "X3096", moreids$Sample.Id)
moreids$Sample.Id <- gsub("-", ".", moreids$Sample.Id)

###plot scores
source("../../bin/wgsAnalysis.R")
source("../../analysis/2016-09-20/geneSampleMatrix.R")

plotPathwayAcrossScores<-function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName){
  ##figure out how many samples we have in common
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to check for mutation correlates'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes[-1]
  if(TRUE %in% index) {
    pathwayMutated <- apply(patient.sample.muts[index, overlap, drop = FALSE], 2, any)
    ##add patient to data frame
    df<-data.frame(pathwayMutated=pathwayMutated ,score=patient.sample.vars[overlap],patient=sapply(overlap,function(x) paste(unlist(strsplit(x,split='tumor'))[1],collapse=' ')))
    print(df)
    pval<-'1.0'
    ##calculate p-value
    try(pval<-format(t.test(score~pathwayMutated,data=df)$p.value,digits=4), silent = TRUE)
    ##make color by patient and shape by mutation
    p<-ggplot(data=df,aes(x=pathwayMutated,y=score))+geom_boxplot(aes(x=pathwayMutated,y=score), outlier.shape = NA)+geom_jitter(aes(x=pathwayMutated,y=score,color=patient))+
      ggtitle(paste(paste(pathwayGenes[1],collapse="_"),'mutation status by\n',testName,'scores'))
    ggsave(paste(paste(pathwayGenes[1],collapse="_"),'mutationsBy',testName,'score.png',sep='_'))
    file.name <- paste(pathwayGenes[1],'mutationsBy',testName,'score.png',sep='_')
    #synStore(File(file.name, parentId='syn7874685'), executed=this.file)
    return(df)
  } else {
    print(paste(pathwayGenes[1], " not in seq data"))
  }
}

calculatePvalues <- function(patient.sample.vars,patient.sample.muts,pathwayGenes,testName) {
  overlap<-intersect(names(patient.sample.vars),colnames(patient.sample.muts))
  print(paste('We have',length(overlap),'samples to evaluate significant difference'))
  index <- row.names(patient.sample.muts) %in% pathwayGenes[-1]
  pathwayMutated <- apply(patient.sample.muts[index, overlap, drop = FALSE], 2, any)
  df<-data.frame(pathwayMutated = pathwayMutated, score = patient.sample.vars[overlap])
  pval<-1
  try((pval <- t.test(df$score ~ df$pathwayMutated)$p.value), silent = TRUE)
  print(pval)
}

som.vars<-somaticGeneSampleMatrix()
germ.vars<-germlineGeneSampleMatrix()

colnames(som.vars) <- gsub(" ", "", colnames(som.vars))
colnames(germ.vars) <- gsub(" ", "", colnames(germ.vars))

##germline variants that correlate with somatic burden....
patients<-sapply(colnames(germ.vars),function(x) unlist(strsplit(x,split='tumor'))[1])

patient.vars<-sapply(unique(patients),function(x){
  apply(germ.vars[,which(patients==x)],1,function(y) any(y))
})

gl.vars.by.sample<-patient.vars[,patients]
colnames(gl.vars.by.sample)<-colnames(germ.vars)

pathway.list <- data.frame(c("CDC27_CREBBP", "CREBBP", "CDC27"), c("CREBBP","CREBBP", NA), c("CDC27","CDC27", NA))

pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

for(x in pathway.list[1]){
for(y in rownames(hallmark.ssGSEA)){
    hall<-as.numeric(hallmark.ssGSEA[y,])
    names(hall)<-colnames(hallmark.ssGSEA)
    ##########
    plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y)))
    names(pvalhallmarkGermline) <- c("Pathway", "p.value")
    pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
    #########
    plotPathwayAcrossScores(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)))
    names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
    pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
    
    pvalhallmarkGermline.df$BH <- p.adjust(pvalhallmarkGermline.df$p.value, method = "BH")
    pvalhallmarkSomatic.df$BH <- p.adjust(pvalhallmarkSomatic.df$p.value, method = "BH")
}}

write.table(pvalhallmarkSomatic.df, "CDC27_CREBBP_pvalhallmarkSomatic.txt")

pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

for(x in pathway.list[2]){
  for(y in rownames(hallmark.ssGSEA)){
    hall<-as.numeric(hallmark.ssGSEA[y,])
    names(hall)<-colnames(hallmark.ssGSEA)
    ##########
    plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y)))
    names(pvalhallmarkGermline) <- c("Pathway", "p.value")
    pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
    #########
    plotPathwayAcrossScores(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)))
    names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
    pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
    
    pvalhallmarkGermline.df$BH <- p.adjust(pvalhallmarkGermline.df$p.value, method = "BH")
    pvalhallmarkSomatic.df$BH <- p.adjust(pvalhallmarkSomatic.df$p.value, method = "BH")
  }}

write.table(pvalhallmarkSomatic.df, "CREBBP_pvalhallmarkSomatic.txt")

pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

for(x in pathway.list[3]){
  for(y in rownames(hallmark.ssGSEA)){
    hall<-as.numeric(hallmark.ssGSEA[y,])
    names(hall)<-colnames(hallmark.ssGSEA)
    ##########
    plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('HallmarkGermline',sep='',y)))
    names(pvalhallmarkGermline) <- c("Pathway", "p.value")
    pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
    #########
    plotPathwayAcrossScores(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('HallmarkSomatic',sep='',y)))
    names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
    pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
    
    pvalhallmarkGermline.df$BH <- p.adjust(pvalhallmarkGermline.df$p.value, method = "BH")
    pvalhallmarkSomatic.df$BH <- p.adjust(pvalhallmarkSomatic.df$p.value, method = "BH")
  }}

write.table(pvalhallmarkSomatic.df, "CDC27_pvalhallmarkSomatic.txt")


cib.scores<-read.table(synGet('syn5809355')@filePath,sep=',',header=T)
cib.scores$Input.Sample<-as.character(cib.scores$Input.Sample)
cib.scores[1,1] <- "3096-PBK-0001"
colnames(cib.scores)[1] <- "Sample.Id"

ids<-read.table(synGet("syn6156140")@filePath, sep = ",", header = T)

colnames(ids)[4] <- "RNASeq (Cufflinks)"
moreids<-synTableQuery("SELECT * FROM syn5556216")@values %>% filter(!is.na(RNASeq)) %>% 
  left_join(ids) %>% 
  dplyr::select(sampleIdentifier, Sample.Id) %>% 
  arrange(Sample.Id)

cib.scores<-cib.scores %>% left_join(moreids)
rownames(cib.scores) <- cib.scores$sampleIdentifier

cib.scores$patient <- t(as.data.frame(strsplit(cib.scores$sampleIdentifier, split = "tumor"))[1,])

celltype <- c("Mast.cells.resting", "Macrophages.M2")


pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

for(x in pathway.list[1]){
  for(y in celltype){
    hall<-as.numeric(cib.scores[,y])
    names(hall)<-cib.scores$sampleIdentifier
    print(hall)
    ##########
    plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('CIBERSORTGermline',sep='',y))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('CIBERSORTGermline',sep='',y)))
    names(pvalhallmarkGermline) <- c("Pathway", "p.value")
    pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
    #########
    plotPathwayAcrossScores(hall,som.vars,x,paste('CIBERSORTSomatic',sep='',y)) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('CIBERSORTSomatic',sep='',y)))
    names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
    pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
    
    pvalhallmarkGermline.df$BH <- p.adjust(pvalhallmarkGermline.df$p.value, method = "BH")
    pvalhallmarkSomatic.df$BH <- p.adjust(pvalhallmarkSomatic.df$p.value, method = "BH")
  }}

write.table(pvalhallmarkSomatic.df, "CREBBP_CDC27_pvalCIBERSORTSomatic.txt")

pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

for(x in pathway.list[2]){
  for(y in celltype){
    hall<-as.numeric(cib.scores[,y])
    names(hall)<-cib.scores$sampleIdentifier
    print(hall)
    ##########
    plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('CIBERSORTGermline',sep='',y))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('CIBERSORTGermline',sep='',y)))
    names(pvalhallmarkGermline) <- c("Pathway", "p.value")
    pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
    #########
    plotPathwayAcrossScores(hall,som.vars,x,paste('CIBERSORTSomatic',sep='',y)) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('CIBERSORTSomatic',sep='',y)))
    names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
    pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
    
    pvalhallmarkGermline.df$BH <- p.adjust(pvalhallmarkGermline.df$p.value, method = "BH")
    pvalhallmarkSomatic.df$BH <- p.adjust(pvalhallmarkSomatic.df$p.value, method = "BH")
  }}

write.table(pvalhallmarkSomatic.df, "CREBBP_pvalCIBERSORTSomatic.txt")

pvalhallmarkGermline.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)
pvalhallmarkSomatic.df <- data.frame("Pathway"=c("placeholder"), "p.value"=c(1), stringsAsFactors = FALSE)

for(x in pathway.list[3]){
  for(y in celltype){
    hall<-as.numeric(cib.scores[,y])
    names(hall)<-cib.scores$sampleIdentifier
    print(hall)
    ##########
    plotPathwayAcrossScores(hall,gl.vars.by.sample,x,paste('CIBERSORTGermline',sep='',y))
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkGermline <- c(paste(pathwayname," ",y), calculatePvalues(hall,gl.vars.by.sample,x,paste('CIBERSORTGermline',sep='',y)))
    names(pvalhallmarkGermline) <- c("Pathway", "p.value")
    pvalhallmarkGermline.df <- rbind(pvalhallmarkGermline.df, pvalhallmarkGermline)
    #########
    plotPathwayAcrossScores(hall,som.vars,x,paste('CIBERSORTSomatic',sep='',y)) 
    pathwayname <- capture.output(print(x[1], max.levels = 0))
    pathwayname <- sub(pattern = "[1]", pathwayname, replacement = "", fixed = TRUE)
    pvalhallmarkSomatic <- c(paste(pathwayname," ",y), calculatePvalues(hall,som.vars,x,paste('CIBERSORTSomatic',sep='',y)))
    names(pvalhallmarkSomatic) <- c("Pathway", "p.value")
    pvalhallmarkSomatic.df <- rbind(pvalhallmarkSomatic.df, pvalhallmarkSomatic)
    
    pvalhallmarkGermline.df$BH <- p.adjust(pvalhallmarkGermline.df$p.value, method = "BH")
    pvalhallmarkSomatic.df$BH <- p.adjust(pvalhallmarkSomatic.df$p.value, method = "BH")
  }}

write.table(pvalhallmarkSomatic.df, "CDC27_pvalCIBERSORTSomatic.txt")

