#plot GSVA scores
library(synapseClient)
library(pheatmap)
source("../../bin/wgsAnalysis.R")
som.vars<-somaticGeneSampleMatrix()
this.script='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-12-19/rePlotGsva.R'

gsva.hallmarks<-read.table(synGet('syn7862153')@filePath)
#gsea.hallmarks<-read.table(synGet('syn7862155')@filePath)

colnames(gsva.hallmarks)<-sapply(colnames(gsva.hallmarks),function(x) tolower(gsub('.',' ',x,fixed=T)))

gsva.hallmarks<-gsva.hallmarks[,grep('patient',colnames(gsva.hallmarks))]
overlap<-intersect(colnames(gsva.hallmarks),colnames(som.vars))
log.df<-data.frame(CDC27=som.vars['CDC27',overlap],CREBBP=som.vars['CREBBP',overlap])
ord<-order(rowSums(log.df*1))
#gh<-gsva.hallmarks[,ord]
df<-data.frame(apply(log.df,2,as.factor))

has.mut.data<-data.frame(HasMutationData=as.factor(sapply(colnames(gsva.hallmarks),function(x) x%in%rownames(log.df))))

pheatmap(gsva.hallmarks[,overlap][,ord],cluster_cols =F,cellheight=10,cellwidth = 10,annotation_col = df,file='gsvaWithCrebbpCDC27.png')
pheatmap(gsva.hallmarks,cluster_cols =T,cellheight=10,cellwidth = 10,annotation_col = has.mut.data,file='gsvaAll.png')

synStore(File('gsvaAll.png',parentId='syn7818711'),used=list(list(entity='syn7862153')),executed=list(list(url=this.script)))
synStore(File('gsvaWithCrebbpCDC27.png',parentId='syn7818711'),used=list(list(entity='syn7862153')),executed=list(list(url=this.script)))