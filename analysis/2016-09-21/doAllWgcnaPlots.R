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
  #  cl.enrich=getEnrichment(cl.clust$expr,cl.clust$tomStatic,cl.clust$TOMprefix,ntop=20)
    cl.eigen=evalEigenModules(cl.clust$expr,colorh1=cl.clust$tomStatic,pids=patient_tumor_number_rna(rownames(cl.clust$expr)),prefix=cl.clust$TOMprefix)

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


createMutMatrix<-function(mut.matrix,genes,patients){
    genes<-intersect(genes,rownames(mut.matrix))
    mmat<-mut.matrix[genes,]

    if(length(genes)==1){
        mmat<-t(mmat)
    }
    rownames(mmat)<-genes

    num.vars<-sapply(tolower(patients),function(p){
        if(p%in%colnames(mmat)){
            nvec=rep(0,nrow(mmat))
            nvec[which(mmat[,p])]<-1

        }else{
            nvec=rep(2,nrow(mmat))

        }
        return(nvec)})
    if(length(genes)==1)
        num.vars<-t(num.vars)
    rownames(num.vars)<-genes
    colnames(num.vars)<-patients
    return(num.vars)

}

modules<-colnames(cl.eigen)
#pdf('TOMEigenGenePlots.pdf')
for(m in modules){
                                        # sizeGrWindow(8,7)
    gl.genes<-unlist(strsplit(enriched.genes[m,'germline'],split=','))
    som.genes<-unlist(strsplit(enriched.genes[m,'somatic'],split=','))

    #now get the patients for which those mutations exists
    which.module<-gsub("ME","",m)
    
    pdf(paste(which.module,'moduleInGeneExpressionAndMuts.pdf',sep=''),height=18)
    
    ME=cl.eigen[,m]
    par(mfrow=c(ifelse(length(som.genes)==0,3,4),1), mar=c(0.3, 5.5, 8, 2))
    plotMat(t(scale(datExpr[order(ME),colorh1==which.module ])),
            nrgcols=30,rlabels=F,rcols=which.module,clabels=rownames(datExpr)[order(ME)],#sapply(rownames(datExpr),function(x) gsub(' ','_',x)),
            main=which.module, cex.main=2)
    par(mar=c(5, 4.2, 0, 0.7))
    barplot(ME[order(ME)], col=which.module, main="", cex.main=2,
            ylab="eigengene expression",xlab="dermalNF sample")

    ##add heatmap for mutations
    if(length(som.genes)>0){
        som.mat<-createMutMatrix(som.mut.matrix,som.genes,rownames(cl.eigen))[,order(ME)]
        if(length(som.genes)==1)
            som.mat<-t(som.mat)
        rownames(som.mat)<-som.genes
        par( mar=c(0.2, 5.5, 0.2, 2))
        image(1:ncol(som.mat),1:nrow(som.mat),z=t(som.mat),col=c('red','blue','grey'),axes=F,xlab='',ylab='')
       # axis(1, at = 1:ncol(som.mat), labels = colnames(som.mat), las = 2, cex.axis = 0.6, 
      #       col.axis = 1)
        axis(2, at = 1:length(som.genes), labels = som.genes, las = 2, cex.axis = 0.6, 
             col.axis = 1)
    }
    if(length(gl.genes)>0){
        germ.mat<-createMutMatrix(germ.mut.matrix,gl.genes,rownames(cl.eigen))[,order(ME)]
        par( mar=c(8, 5.5, 0.2, 2))
        
        image(1:ncol(germ.mat),1:nrow(germ.mat),z=t(germ.mat),col=c('red','blue','grey'),axes=F,xlab='',ylab='')
        axis(1, at = 1:ncol(germ.mat), labels = colnames(germ.mat), las = 2, cex.axis = 0.6, 
             col.axis = 1)
        axis(2, at = 1:length(gl.genes), labels = gl.genes, las = 2, cex.axis = 0.1, 
             col.axis = 1)     
        
    }
    dev.off()
}
#dev.off()
all.files=list.files('.')
pdf.files=all.files[grep('pdf',all.files)]



#clfiles=c('cuffLinksWGCNAClustering.pdf','WGCNA_cuffLinksTOMClusterAssignment.tsv','cuffLinksTOMtop10GOTermsPermodule.csv')
for(f in pdf.files){
#f='TOMEigenGenePlots.pdf'


synStore(File(f,parentId='syn5669860'),
              used=c(exlist,list(list(entity='syn5579598',wasExecuted=F))))

}

                                        #store on synapse
