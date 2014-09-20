# Comparisons.R
# Copyright 2014 Matt Ralston
# Updated 9/11/2014
#
# This R script contains code to produce tables of pairwise comparisons
# from differential expression

# 2 p-values 0.01, 0.05
# 2 Log2 fold changes 1.5 2.0
# 2 stress levels
# 4 time points
# BuOH table
#
# a = alpha value
# b = log2 fold change
library('gtools')

#       D 3    S T U F F


times<-c("15","75","150","270")
stresses<-c("BA","BuOH","NS")
texes<-c("Normal","TEX")

time1=list()
time2=list()

buohtime=list()
batime=list()
# round to 'r' digits
r = 5





# This compares all possible sets of combinations, printing the results to an svg

for (t1 in times) {
    for (s1 in stresses) {
        # Control conditions
        for (t2 in times) {
            for (s2 in stresses) {
                # Multi factor comparison
                if (s1 != s2 || t1 != t2) {
                    x<-as.data.frame(results(deseq,contrast=list(paste('time',t1,'.cond',s1,sep=''),paste('time',t2,'.cond',s2,sep=''))))
                    l<-length(rownames(x))
                    cndtn<-rep(paste(s1,'.',t1,sep=''),l)
                    cntrl<-rep(paste(s2,'.',t2,sep=''),l)
                    x$padj[is.na(x$padj)]<-0.9999999
                    x<-data.frame(id=rownames(x),Logexp=round(log10(x$baseMean+1),digits=r),foldchange=round(x$log2FoldChange,digits=r),pval=round(pval.trans(x$padj),digits=r),condition=cndtn,control=cntrl)
                    write.table(x,file=paste("summary/data/",s1,'.',t1,'-',s2,'.',t2,".csv",sep=''),quote=F,sep=",",row.names=F,col.names=T)
                }
                # Stress comparison
                if (s1 != s2) {
                    x<-as.data.frame(results(deseq,contrast=list(paste('cond',s1,sep=''),paste('cond',s2,sep=''))))
                    l<-length(rownames(x))
                    cndtn<-rep(s1,l)
                    cntrl<-rep(s2,l)
                    x$padj[is.na(x$padj)]<-0.9999999           
                    x<-data.frame(id=rownames(x),Logexp=round(log10(x$baseMean+1),digits=r),foldchange=round(x$log2FoldChange,digits=r),pval=round(pval.trans(x$padj),digits=r),condition=cndtn,control=cntrl)
                    write.table(x,file=paste("summary/data/",s1,'.',s2,'.csv',sep=''),quote=F,sep=",",row.names=F,col.names=T)
                }
            }
            # Time comparison
            if (t1 != t2) {
                x<-as.data.frame(results(deseq,contrast=list(paste('time',t1,sep=''),paste('time',t2,sep=''))))
                l<-length(rownames(x))
                cndtn<-rep(t1,l)
                cntrl<-rep(t2,l)
                x$padj[is.na(x$padj)]<-0.9999999
                x<-data.frame(id=rownames(x),Logexp=round(log10(x$baseMean+1),digits=r),foldchange=round(x$log2FoldChange,digits=r),pval=round(pval.trans(x$padj),digits=r),condition=cndtn,control=cntrl)
                write.table(x,file=paste("summary/data/",t1,'.',t2,'.csv',sep=''),quote=F,sep=",",row.names=F,col.names=T)
            }
        }
    }
}












