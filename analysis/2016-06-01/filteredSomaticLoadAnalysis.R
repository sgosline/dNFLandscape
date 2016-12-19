##compare TCGA mutational analysis to dermal somatic load

#source("../../bin/dermalNFData.R")

require(dplyr)
require(ggplot2)
source("../../bin/wgsAnalysis.R")

#first do somatic load
##now get the TCGA data
source("../../bin/TcgaMutationalData.R")
all.genes<-read.table('../../data/HugoGIDsToEntrez_DAVID.txt',header=T,as.is=T,sep='\t',quote='"')[,1]

non.silent.tcga<-subset(combinedMaf,Variant_Classification!='Silent')
dft <- non.silent.tcga%>%group_by(Patient)%>%summarize(MutatedGenes=n_distinct(Hugo_Symbol))
dft$Disease<-non.silent.tcga$tumor_type[match(dft$Patient,non.silent.tcga$Patient)]

effects=c("HIGH","HIGH_MODERATE_LOW","HIGH_MODERATE")
soms=c("LikelyLOH_StrongLOH_LikelySomatic_StrongSomatic")
for(e in effects){
  for(s in soms){
    tab<-subset(expr.gene.muts1,Effect%in%unlist(strsplit(e,split='_')))
    tab<-subset(tab,Status%in%unlist(strsplit(s,split='_')))
    scounts<-tab%>%group_by(Sample)%>%summarize(nMuts=n_distinct(Gene))
    
  }
}




dft<-rbind(dft,data.frame(Disease=rep('dermalNF',nrow(som.counts)),Patient=som.counts$Sample_ID,
                          MutatedGenes=som.counts$MutatedGenes))

require(ggplot2)
p<-ggplot(dft)+geom_boxplot(aes(x=Disease,y=MutatedGenes))+scale_y_log10()+theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)

png('somaticBurdenPerGene.png')
print(p)
dev.off()

dfa <- non.silent.tcga%>%group_by(Tumor_Sample_Barcode)%>%summarize(DistinctMutations=n())
dfa$Disease<-non.silent.tcga$tumor_type[match(dfa$Tumor_Sample_Barcode,non.silent.tcga$Tumor_Sample_Barcode)]

dfa<-rbind(dfa,data.frame(Disease=rep('dermalNF',nrow(som.counts)),Tumor_Sample_Barcode=som.counts$Sample_ID,
                          DistinctMutations=som.counts$DistinctMutations))

require(ggplot2)
p<-ggplot(dfa)+geom_boxplot(aes(x=Disease,y=DistinctMutations))+scale_y_log10()+theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)


##what if we filtered for express genes
rna.mat<-rna_count_matrix(doLogNorm=FALSE,minCount=2)
expr.non.silent=subset(non.silent,Hugo_Symbol%in%rownames(rna.mat))

expr.som.counts=subset(expr.non.silent,Mutation_Status=='Somatic') %>% group_by(Sample_ID) %>% summarize(MutatedGenes=n_distinct(Hugo_Symbol),DistinctMutations=n())
edf=tidyr::gather(expr.som.counts,'SomaticEvent','Count',2:3)
p<-ggplot(edf)+geom_bar(aes(x=Sample_ID,y=Count,fill=SomaticEvent),stat='identity',position='dodge')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

edft <- non.silent.tcga%>%group_by(Patient)%>%summarize(MutatedGenes=n_distinct(Hugo_Symbol))
edft$Disease<-non.silent.tcga$tumor_type[match(edft$Patient,non.silent.tcga$Patient)]

edft<-rbind(edft,data.frame(Disease=rep('dermalNF',nrow(expr.som.counts)),Patient=expr.som.counts$Sample_ID,
                          MutatedGenes=expr.som.counts$MutatedGenes))
p<-ggplot(edft)+geom_boxplot(aes(x=Disease,y=MutatedGenes))+scale_y_log10()+theme(axis.text.x=element_text(angle = -90, hjust = 0))
print(p)
png('somaticBurdenPereExpressedGene.png')
print(p)
dev.off()

