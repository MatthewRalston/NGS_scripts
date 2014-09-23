#!/usr/bin/env Rscript
library('cluster')
options(warn=-1)
args<- commandArgs(trailingOnly=TRUE)
clusters<-read.table(file=args[1],sep=",",col.names=c("id","cluster"))
data<-read.table(file=args[2],sep=",")
sils<-silhouette(clusters$cluster,dist(data[clusters$id,]))
cat(sils[,3])


