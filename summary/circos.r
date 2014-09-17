# Circos.R
# Copyright 2014 Matt Ralston
# Updated 9/4/2014
#
# This R script contains code to generate files tables from DESeq2
# Differential expression analyses for processing and eventual input into circos

times<-c("15","75","150","270")
stresses<-c("BA","BuOH","NS")
texes<-c("Normal","TEX")

# Time only
for (i in 2:4) {
    #  Type 1: (ti vs t1)
    #  Type 2: (ti vs t[i-1])
    t1<-results(deseq,contrast=list(paste('time',times[i],sep=''),paste('time',times[1],sep='')))
    t2<-results(deseq,contrast=list(paste('time',times[i],sep=''),paste('time',times[i-1],sep='')))
    t1<-data.frame(GeneID=rownames(t1),Log2FC=t1$log2FoldChange,pValue=t1$padj)
    t2<-data.frame(GeneID=rownames(t2),Log2FC=t2$log2FoldChange,pValue=t2$padj)
    write.table(t1,file=paste("circos/data/foldchanges/time/typea/",times[i],"v15.pre",sep=''),sep="\t",col.names=FALSE, row.names=FALSE,quote=FALSE)
    write.table(t2,file=paste("circos/data/foldchanges/time/typeb/",times[i],'v',times[i-1],".pre",sep=''),sep="\t",col.names=FALSE, row.names=FALSE,quote=FALSE)
}

for (j in 1:4) {
    butyrate<-results(deseq,contrast=list(paste('time',times[j],'.cond',stresses[1],sep=''),paste('time',times[i],'.cond',stresses[3],sep='')))
    butanol<-results(deseq,contrast=list(paste('time',times[j],'.cond',stresses[2],sep=''),paste('time',times[i],'.cond',stresses[3],sep='')))
    butyrate<-data.frame(GeneID=rownames(butyrate),Log2FC=butyrate$log2FoldChange,pValue=butyrate$padj)
    butanol<-data.frame(GeneID=rownames(butanol),Log2FC=butanol$log2FoldChange,pValue=butanol$padj)
    write.table(butyrate,file=paste("circos/data/foldchanges/stress/butyrate",times[j],".pre",sep=''),sep="\t",col.names=FALSE, row.names=FALSE,quote=FALSE)
    write.table(butanol,file=paste("circos/data/foldchanges/stress/butanol",times[j],".pre",sep=''),sep="\t",col.names=FALSE, row.names=FALSE,quote=FALSE)
}

ba<-results(deseq,contrast=list('condBA','condNS'))
buoh<-results(deseq,contrast=list('condBuOH','condNS'))
ba<-data.frame(GeneID=rownames(ba),Log2FC=ba$log2FoldChange,pValue=ba$padj)
buoh<-data.frame(GeneID=rownames(buoh),Log2FC=buoh$log2FoldChange,pValue=buoh$padj)
write.table(ba,file="circos/data/foldchanges/stress/ba.pre",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(buoh,file="circos/data/foldchanges/stress/buoh.pre",sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
