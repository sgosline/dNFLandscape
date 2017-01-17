##will take a long time to run
source('SNP_conversion.R')
library(MatrixEQTL)
library(synapseClient)
library(tidyr)
library(dplyr)
used <- c("syn5492805", "syn5556216", "syn5004874", "syn7342635", "syn5297573", "syn5004912", "syn8030595")

useModel = modelLINEAR
SNP_file_name = "SNPsinRNA.txt"
expression_file_name = "RNAinSNPs.txt"
covariates_file_name = "covariates.txt"
output_file_name_cis = "output_cis.txt"
output_file_name_tra = "output_tra.txt"
pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-2
errorCovariance = numeric()
cisDist = 1e6

snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
  snps$LoadFile( SNP_file_name )

gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
  gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  cvrt$fileSliceSize = 2;      # read file in pieces of 2,000 rows
  cvrt$LoadFile( covariates_file_name );

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = 0, ##allows running in cis-mode only - can run trans-mode too by turning on next line and turning this off
  #pvOutputThreshold.tra = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

cis <- me$cis$eqtls

plot(me, pch = 16, cex = 0.7)

me2 = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = 0, ##allows running in cis-mode only - can run trans-mode too by turning on next line and turning this off
  #pvOutputThreshold.tra = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(me2, col='grey')


##reference snp names
rsid.index<-fread(file="InfiniumOmni2-5-8v1-3_A1_b144_rsids.txt", head = TRUE, sep="\t", stringsAsFactors = FALSE)
rsid.index<-filter(rsid.index, RsID!=".")
colnames(rsid.index)[1] <- 'snps'

cis<- left_join(cis,rsid.index,by = "snps")
cis<-filter(noX_Y.df, FDR<0.05)

siggene<-count(cis, gene)
sigsnps<-count(cis, snps)

write.table(siggene, file = "sig_genes.txt")
write.table(sigsnps, file = "sig_snps.txt")

##use snpspos_noXY for cis eQTL analysis without X or Y chromosome position data, so these chromosomes will not be tested
output_file_name_cis = 'output_cis_noXY.txt'

noX_Y = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = 0, ##allows running in cis-mode only - can run trans-mode too by turning on next line and turning this off
  #pvOutputThreshold.tra = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos_noXY, 
  genepos = genepos_noXY,
  cisDist = cisDist,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(noX_Y, col='grey')
noX_Y.df<- noX_Y$cis$eqtls
noX_Y.df<- left_join(noX_Y.df,rsid.index,by = "snps")
noX_Y.df<-filter(noX_Y.df, FDR<0.05)

siggene<-count(noX_Y.df, gene)
sigsnps<-count(noX_Y.df, snps)

write.table(siggene, file = "sig_genes_noXY.txt")
write.table(sigsnps, file = "sig_snps_noXY.txt")

noX_Y_2 = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = 0, ##allows running in cis-mode only - can run trans-mode too by turning on next line and turning this off
  #pvOutputThreshold.tra = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos_noXY, 
  genepos = genepos_noXY,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

plot(noX_Y_2, pch = 16, cex = 0.7)

##plot snp measurement vs gene expression
this.file='https://raw.githubusercontent.com/allaway/dNFLandscape/master/analysis/2016-12-28/matrix_eQTL.R'
library(ggplot2)

PlotSigSnps<-function(df, nrows, nameofdata){
  snps<-as.character(df$snps)
  genes<-as.character(df$gene)
  pvalue<-df$FDR
  rsID<-as.character(df$RsID)
for(x in c(1:nrows)){
  chosen.snp <- snps[x]
  chosen.gene <- genes[x]
  FDR <- signif(pvalue[x], 3)
  this.rsid <- rsID[x]
  snpdat <-as.numeric(filter(snp2, id==chosen.snp)[,-1])
  genedat <-as.numeric(filter(expression2, id==chosen.gene)[,-1])
  qplot(x = snpdat, y = genedat, color = snpdat, size = genedat, xlim = c(-0.5, 2.5), ylab = "log(2)FPKM", xlab = "SNP_Measurement",  
        main = paste(chosen.snp,"_",chosen.gene,"_expression_level_FDR = ",FDR,"_rsID=",this.rsid))
  ggsave(filename = paste(chosen.snp,"_",this.rsid,"_",chosen.gene,"_expr_",nameofdata,".png",  sep = ""))
  synStore(File(paste(chosen.snp,"_",this.rsid,"_",chosen.gene,"_expression_level_",nameofdata,".png",  sep = ""), parentId='syn8030585'), executed=this.file, used = used)
  }
}


X_Y_incl <- cis
PlotSigSnps(noX_Y.df, nrow(noX_Y.df), "noX_Y")
PlotSigSnps(X_Y_incl, nrow(X_Y_incl), "X_Y_incl")

tumorsize <- select(sample.info, Patient, TumorNumber, Length_in_mm)
tumorsize <- unite(tumorsize, "PatientTumorNumber", Patient, TumorNumber, sep = "-", remove = TRUE)

sizes <- filter(tumorsize, PatientTumorNumber %in% colnames(expression2))

PlotSigSnpsAndTumorSize<-function(df, nrows, nameofdata){
  snps<-as.character(df$snps)
  genes<-as.character(df$gene)
  pvalue<-df$FDR
  rsID<-as.character(df$RsID)
  for(x in c(1:nrows)){
    chosen.snp <- snps[x]
    chosen.gene <- genes[x]
    FDR <- signif(pvalue[x], 3)
    snpdat <- as.numeric(filter(snp2, id==chosen.snp)[,-1])
    genedat <- as.numeric(filter(expression2, id==chosen.gene)[,-1])
    samps <- colnames(expression2)[-1]
    sizes <- filter(tumorsize, PatientTumorNumber %in% samps)
    size_in_mm <- as.numeric(sizes$Length_in_mm)
    qplot(x = snpdat, y = genedat, color = snpdat, size = size_in_mm, xlim = c(-0.5, 2.5), ylab = "log(2)FPKM", xlab = "SNP_Measurement",  
          main = paste(chosen.snp,"_",chosen.gene,"_expression_level_FDR = ",FDR,"_rsID=",this.rsid))
    ggsave(filename = paste(chosen.snp,"_",this.rsid,"_",chosen.gene,"_expr_",nameofdata,"tumor_size_plot.png", sep = ""))
    synStore(File(paste(chosen.snp,"_",this.rsid,"_",chosen.gene,"_expression_level_",nameofdata,"tumor_size_plot.png", sep = ""), parentId='syn8030585'), executed=this.file, used = used)
  }
}

PlotSigSnpsAndTumorSize(noX_Y.df, nrow(noX_Y.df), "noX_Y")
PlotSigSnpsAndTumorSize(X_Y_incl, nrow(X_Y_incl), "X_Y_incl")
