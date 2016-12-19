##compare TCGA mutational analysis to dermal somatic load

#source("../../bin/dermalNFData.R")

require(dplyr)
require(ggplot2)
source("../../bin/wgsAnalysis.R")

#first do somatic load
##now get the TCGA data
source("../../bin/TcgaMutationalData.R")
all.genes<-read.table('../../../dermalNF/data/HugoGIDsToEntrez_DAVID.txt',header=T,as.is=T,sep='\t',quote='"')[,1]

non.silent.tcga<-subset(combinedMaf,Variant_Classification!='Silent')
dft <- non.silent.tcga%>%group_by(Patient)%>%summarize(MutatedGenes=n_distinct(Hugo_Symbol))
dft$Disease<-non.silent.tcga$tumor_type[match(dft$Patient,non.silent.tcga$Patient)]

effects=c("HIGH","HIGH_MODERATE_LOW","HIGH_MODERATE")
soms=c("LikelyLOH_StrongLOH_LikelySomatic_StrongSomatic","LikelySomatic_StrongSomatic","StrongLOH_StrongSomatic")
require(ggplot2)

for(e in effects){
      for(s in soms){
              tab<-subset(expr.gene.muts1,Effect%in%unlist(strsplit(e,split='_')))
                  tab<-subset(tab,Status%in%unlist(strsplit(s,split='_')))
                      scounts<-tab%>%group_by(Sample)%>%summarize(nMuts=n_distinct(Gene))
                          dft2<-rbind(dft,data.frame(Disease=rep('dermalNF',nrow(scounts)),Patient=scounts$Sample,
                                                        MutatedGenes=scounts$nMuts))
                          
                            p<-ggplot(dft2)+geom_boxplot(aes(x=reorder(Disease,MutatedGenes,FUN=median),y=MutatedGenes))+scale_y_log10()+theme(axis.text.x=element_text(angle = -90, hjust = 0))
                            print(p)
                                      
                                          png(paste('somaticBurden',e,s,'PerGene.png',sep=''))
                                              print(p)
                                                  dev.off()
                                                    }
}

