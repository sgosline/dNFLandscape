source("../../bin/clusterRNASeqData.R")
source("../2016-09-20/geneSampleMatrix.R")

source("../../dermalNF/bin/dermalNFData.R")

exlist=list(list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/analysis/2016-09-21/doAllWgcnaPlots.R',wasExecuted=TRUE),
            list(url='https://raw.githubusercontent.com/sgosline/dNFLandscape/master/bin/clusterRNASeqData.R',wasExecuted=TRUE))

#for(signed in c(FALSE,TRUE)){
#  for(sd in c(FALSE,TRUE)){
signed=TRUE
sd=FALSE
#   fc.clust=clusterData(t(fc.matrix),'featureCounts',signed,sd,topGenes=3000)
  #  fc.enrich=getEnrichment(fc.clust$expr,fc.clust$tomStatic,fc.clust$TOMprefix)
#    fc.eigen=evalEigenModules(fc.clust$expr,colorh1=fc.clust$tomStatic,pids=fc.pids,prefix=fc.clust$TOMprefix)
#    write.table(fc.eigen,paste('featureCounts',ifelse(sd,'sd','conn'),'filtered',ifelse(signed,'signed','unsigned'),'clusterEigenGenes.tab',sep='_'))

    ##then get clusters for cufflinks
    cl.clust=clusterData(t(cl.matrix),'cuffLinks',signed,sd,topGenes=3000)
    ##now plot eigen genes for each
    cl.enrich=getEnrichment(cl.clust$expr,cl.clust$tomStatic,cl.clust$TOMprefix,ntop=20)
    cl.eigen=evalEigenModules(cl.clust$expr,colorh1=cl.clust$tomStatic,pids=cl.pids[rownames(cl.clust$expr)],prefix=cl.clust$TOMprefix)

rownames(cl.eigen)<-patient_tumor_number_rna(rownames(cl.eigen))
datExpr<-cl.clust$expr
rownames(datExpr)<-patient_tumor_number_rna(rownames(datExpr))
colorh1<-cl.clust$tomStatic
                                        #write.table(cl.eigen,paste('cuffLinks',ifelse(sd,'sd','conn'),'filtered',ifelse(signed,'signed','unsigned'),'clusterEigenGenes.tab',sep='_'))
#}


                                        #let's get the somatic/germline variants in there too
enriched.genes<-read.table(synGet('syn7256514')@filePath)
##also get mutation matrix
som.mut.matrix<-somaticGeneSampleMatrix()
germ.mut.matrix<-germlineGeneSampleMatrix()

modules<-colnames(cl.eigen)
pdf('TOMEigenGenePlots.pdf')
for(m in modules){
                                        # sizeGrWindow(8,7)
    gl.genes<-unlist(strsplit(enriched.genes[m,'germline'],split=','))
    som.genes<-unlist(strsplit(enriched.genes[m,'somatic'],split=',')))

    #now get the patients for which those mutations exists

    ME=cl.eigen[,m]
    par(mfrow=c(2,1), mar=c(0.3, 5.5, 8, 2))
    which.module<-gsub("ME","",m)
    plotMat(t(scale(datExpr[order(ME),colorh1==which.module ])),
            nrgcols=30,rlabels=F,rcols=which.module,clabels=rownames(datExpr),#sapply(rownames(datExpr),function(x) gsub(' ','_',x)),
            main=which.module, cex.main=2)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME[order(ME)], col=which.module, main="", cex.main=2,
            ylab="eigengene expression",xlab="dermalNF sample")

    ##add heatmap for mutations

}
dev.off()
#all.files=list.files('.')
#pdf.files=all.files[grep('TOM',all.files)]



#clfiles=c('cuffLinksWGCNAClustering.pdf','WGCNA_cuffLinksTOMClusterAssignment.tsv','cuffLinksTOMtop10GOTermsPermodule.csv')
#for(f in pdf.files){
f='TOMEigenGenePlots.pdf'


synStore(File(f,parentId='syn5669860'),
              used=c(exlist,list(list(entity='syn5579598',wasExecuted=F))))

#}

                                        #store on synapse
