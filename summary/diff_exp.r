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



Cairo()
mydir="/home/mrals/Final"
setwd(mydir)
source("summary/functions.r")


#################################################################
#            Data loading and DESeq2 calculations
#################################################################
directory<-"Expression/counts"
sampleFiles<-list.files(paste(mydir,directory,sep="/"),pattern="*.counts.txt",full.names=T)
sampleTime<-rep(c(150,15,270,270,75,75,150,15,270,75),3)
sampleRep<-rep(c(rep("A",6),rep("B",4)),3)
sampleCond<-c(rep("BA",10),rep("BuOH",10),rep("NS",10))
sampleTex<-rep(c(rep("Normal",3),"TEX","Normal","TEX",rep("Normal",4)),3)
sampleNames<-gsub("-","_",unlist(lapply(strsplit(unlist(lapply(strsplit(sampleFiles,"/"),"[[",7)),"\\."),"[[",1)))
samples<-data.frame(sampleName=sampleNames,fileName=sampleFiles,cond=sampleCond,time=sampleTime,replicate=sampleRep,tex=sampleTex)
htseq<-DESeqDataSetFromHTSeqCount(sampleTable=samples,directory='',design= ~ cond + time + replicate + tex)

deseq<-DESeq(htseq)
rdeseq<-rlogTransformation(deseq,blind=TRUE)
rawcounts<-as.data.frame(counts(deseq,normalized=FALSE))
normcounts<-as.data.frame(counts(deseq,normalized=TRUE))
regcounts<-as.data.frame(assay(rdeseq))
#################################################################
#            Exploratory plots
#################################################################
mraw<-melt(rawcounts);mnorm<-melt(normcounts);mreg<-melt(regcounts);
source("summary/correlation.r")

