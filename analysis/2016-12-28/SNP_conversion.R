library(limma)
library(biomaRt)
library(dplyr)
library(tidyr)
library(data.table)
library(synapseClient)

##get gene expression data 
##get dermal NF data 
source("../../dermalNF/bin/dermalNFData.R")
dermals=rna_fpkm_matrix(byIsoform=FALSE)

##step 1 - just combine all
comb=dermals

##step 2, remove values below a particular FPKM, let's say 0.1 - for ALL genes
expr<-1:nrow(comb) #which(apply(comb,1,function(x) all(x>0.2)))
expr<-setdiff(expr,expr[union(grep('MIR',rownames(comb)[expr]),grep("SNO",rownames(comb)[expr]))])

##step 3, normalize
comb.norm=data.frame(voomWithQualityWeights(comb[expr,])$E)
comb.norm$Gene=rownames(comb.norm)

##set up expression data frame
expression <- cbind("id"=comb.norm$Gene, comb.norm[,-34])
colnames(expression) <- gsub("X", "", colnames(expression))
colnames(expression) <- gsub("\\.", "-", colnames(expression))

##rename and find overlapping samples
res<-synQuery("select patientID,tissueID,sampleID from entity where parentId=='syn5492805'")
tres<-synTableQuery("SELECT Patient,RnaID,TumorNumber,'RNASeq (Cufflinks)' FROM syn5556216 where RnaID is not NULL")@values
idx<-match(res$entity.id,tres$`RNASeq (Cufflinks)`)
dres<-res[which(!is.na(idx)),]
tres<-tres[idx[which(!is.na(idx))],]
full.map.rna<-cbind(dres,tres)
full.map.rna<-unite(full.map.rna, PatientTumorNumber, Patient, TumorNumber, sep = "-")

res<-synQuery("select patientID,tissueID,sampleID from entity where parentId=='syn5004874'")
tres<-synTableQuery("SELECT Patient,DnaID,TumorNumber,SNPArray FROM syn5556216 where RnaID is not NULL")@values
idx<-match(res$entity.id,tres$`SNPArray`)
dres<-res[which(!is.na(idx)),]
tres<-tres[idx[which(!is.na(idx))],]
full.map.snp<-cbind(dres,tres)
full.map.snp<-unite(full.map.snp, PatientTumorNumber, Patient, TumorNumber, sep = "-")

full.map<-right_join(full.map.rna, full.map.snp, by = "PatientTumorNumber")

##find overlapping RNA, rename
expression.samples <- intersect(colnames(expression), full.map$entity.sampleID.x)
expression2 <- expression[, expression.samples]
setnames(expression2, old = full.map$entity.sampleID.x, new = full.map$PatientTumorNumber)
expressionid <- rownames(expression2)
expression2 <- cbind(expressionid, expression2)
colnames(expression2)[1] <- c("id")
write.table(expression2, file="RNAinSNPs.txt", sep='\t', row.names = FALSE)

##set up snp datafile
snp <- fread(input = '3096-PBK_FinalReport.txt', head = TRUE, sep = "\t", skip = 9L)
colnames(snp) <- make.names(colnames(snp))
snp <- dplyr::select(snp, SNP.Name, Sample.ID, B.Allele.Freq)
snp <- spread(snp, Sample.ID, B.Allele.Freq)
colnames(snp)[1] <- "id"
snpid<-snp[,1]
snp <- snp[,-1]*2
snp.samples <- intersect(colnames(snp), full.map$entity.sampleID.y)
snp2 <- snp[, snp.samples, with = FALSE]
setnames(snp2,  old = full.map$entity.sampleID.y, new = full.map$PatientTumorNumber)
snp2 <- cbind(snpid, snp2)
write.table(snp2, file = 'SNPsInRNA.txt', sep = "\t", row.names=FALSE)

##set up genepos data frame
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genepos <- getBM(attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"), filters = "hgnc_symbol",
                 values = rownames(comb.norm), mart = mart)   
colnames(genepos) <- c("geneid", "chr", "s1", "s2")
genepos_noXY <- filter(genepos, chr!='X', chr!='XY', chr!='Y')

##set up snpspos data frame 
conversion.df <- fread(input = 'HumanOmni2-5-8-v1-1-C.csv', head = TRUE, skip=7L, sep=",")
conversion.df <- dplyr::select(conversion.df, Name, Chr, MapInfo)
colnames(conversion.df) <- c("snp", "chr", "pos")
chr <- conversion.df$chr
conversion.df$chr<-chr
snpspos<-as.data.frame(conversion.df)
snpspos_noXY <- filter(snpspos, chr!='X', chr!='XY', chr!='Y')

##covariates
##get cutaneous NF sample information
sample.info<-synTableQuery('SELECT * FROM syn5556216')@values
patient.info<-synTableQuery('SELECT * FROM syn7342635')@values
##prep tumor variables for comparison
sample.info<-dplyr::full_join(sample.info,patient.info, by="Patient")
sample.info<-dplyr::filter(sample.info, !is.na(WGS))
sample.info<-dplyr::filter(sample.info, !is.na(TumorNumber))

covariates<-dplyr::select(sample.info, Patient, TumorNumber, Race, Gender, Age)
covariates$Gender[covariates$Gender=='Female'] <- 1
covariates$Gender[covariates$Gender=='Male'] <- 0
covariates$PatientTumorNumber <- paste(covariates$Patient, covariates$TumorNumber, sep="-")
covariates <- covariates[,3:6]
full.map <- dplyr::left_join(full.map, covariates, by="PatientTumorNumber")
covariates <- dplyr::select(full.map, PatientTumorNumber, Gender, Age)
rownames(covariates) <- covariates$PatientTumorNumber
covariates <- t(covariates[,2:3])
id <- rownames(covariates)
covariates <- cbind(id, covariates)
write.table(covariates, file = 'covariates.txt', sep = "\t", row.names=FALSE, quote = FALSE)
