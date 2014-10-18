# Diff_exp.R
# Copyright 2014 Matt Ralston
# Updated 7/25/2014
#
# This R script contains code for summary and diagnostic plots
# for differential expression data (read counts per gene) from
# RNA-seq data.
library('ggplot2')
library('DESeq2')
library('scales')
library('RColorBrewer')
library('reshape2')
library('cairoDevice')
library('gridExtra')
options(max.print=1000)
options(scipen=999)

Cairo()
mydir="/home/mrals/Final"
setwd(mydir)
source("summary/functions.r")


#################################################################
#            Data loading and DESeq2 calculations
#################################################################
directory<-"Expression/counts"
sampleFiles<-list.files(paste(mydir,directory,sep="/"),pattern="*.counts.txt",full.names=T)
sampleTime<-factor(rep(c(150,15,270,270,75,75,150,15,270,75),3))
sampleRep<-rep(c(rep("A",6),rep("B",4)),3)
sampleCond<-c(rep("BA",10),rep("BuOH",10),rep("NS",10))
sampleTex<-rep(c(rep("Normal",3),"TEX","Normal","TEX",rep("Normal",4)),3)
sampleNames<-gsub("-","_",unlist(lapply(strsplit(unlist(lapply(strsplit(sampleFiles,"/"),"[[",7)),"\\."),"[[",1)))
samples<-data.frame(sampleName=sampleNames,fileName=sampleFiles,cond=sampleCond,time=sampleTime,replicate=sampleRep,tex=sampleTex)
htseq<-DESeqDataSetFromHTSeqCount(sampleTable=samples,directory='',design= ~ replicate + tex + time + cond + time:cond)

deseq<-DESeq(htseq)
rdeseq<-rlogTransformation(deseq,blind=TRUE)
rawcounts<-as.data.frame(counts(deseq,normalized=FALSE))
normcounts<-as.data.frame(counts(deseq,normalized=TRUE))
regcounts<-as.data.frame(assay(rdeseq))
write.table(regcounts,file="circos/data/expression.txt",row.names=T,col.names=T,quote=FALSE)
res<-results(deseq)
#################################################################
#            Exploratory plots
#################################################################
mraw<-melt(rawcounts);mnorm<-melt(normcounts);mreg<-melt(regcounts);
# Correlations [raw, normalized, regularized]
source("summary/correlation.r")
# Jitter plots showing normalization
p1<-ggplot(melt(rawcounts[,c(21,27,22,28)]+1),aes(x=variable,y=value))+geom_jitter()+scale_y_log10(breaks=10**(0:4),labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='l')+stat_summary(fun.y=median,fun.ymax=top.quartile,fun.ymin=bottom.quartile,geom='crossbar',colour='red')+ylab('Counts')+theme(axis.title.x=element_blank(),axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+ggtitle("Raw Counts")
p2<-ggplot(melt(normcounts[,c(21,27,22,28)]+1),aes(x=variable,y=value))+geom_jitter()+scale_y_log10(breaks=10**(0:4),labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='l')+stat_summary(fun.y=median,fun.ymax=top.quartile,fun.ymin=bottom.quartile,geom='crossbar',colour='red')+ylab('Normalized Counts')+theme(axis.title.x=element_blank(),axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+ggtitle("Normalized Counts")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave("summary/images/Normalization.png",p)



# MA-plots and lists of differential expression results (time1 time2 stress buohtime batime)
# for buohtime and batime, indexes 1, 2, 3, 4 each refere to timepoints 15, 75, 150, and 270, resp.
source("summary/maplots.r")
# Number of genes in each comparison
source("summary/comparisons.r")
# Variance vs. expression level trend, Residual plot, dispersion
source("summary/regularization.r")
# PCA
source("summary/pca.r")
# Cluster analysis
source("summary/cluster.r")
# GO analysis


# Circos data
# Currently: logit transform
# other ideas:
# angular -> asin(sqrt(x))
# http://www.stata.com/users/njc/topichlp/transint.hlp
source("summary/circos.r")
