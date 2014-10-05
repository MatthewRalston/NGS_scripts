# Cluster.R
# Copyright 2014 Matt Ralston
# Updated 7/25/2014
#
# This is an R script performs clustering analysis

library('gplots')
library('fpc')
library('cluster')
a = 0.05 # alpha value
l = 2.0 # log fold change threshold


n=100 # of genes
#k=floor(sqrt(n/2)) # rule of thumb

# 0.1, 0.5 -> ba=1473  buoh=603
# 0.05, 1 -> ba=577 148
# 0.1 1 -> ba=632 173



#sig<-subset(deres,deres$padj < a)
#n=length(rownames(sig))

#ordered<-as.vector(rownames(deres[order(deres$padj),]))
#select<-ordered[1:n]
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# test
#heatmap.2(as.matrix(normcounts[select,]), col=hmcol, Rowv=FALSE, Colv=FALSE,scale="none",dendrogram="none",trace="none",margin=c(10,6))
# Distance matrix / heatmap of distances f regularized data
#dists<-dist(t(assay(rdeseq)))
#mat<-as.matrix(dists)
#rownames(mat)<-colnames(mat)<-with(colData(deseq),paste(cond,time,replicate,sep=" : "))
#heatmap.2(mat, trace="none",col=rev(hmcol),margin=c(10,10))
# Clustering

# Data as matrix
#x<-as.matrix(regcounts) # whole dataset
#rownames(x)<-rownames(regcounts)
#x<-as.matrix(regcounts[select,]) #subset of select


#   D A T A      P R E P
ba<-2^regcounts[rownames(subset(stress[[2]],subset=stress[[2]]$padj < a & abs(stress[[2]]$log2FoldChange) > l)),]
ba<-(ba-rowMeans(ba))/apply(ba,1,sd)
nba=length(rownames(ba))
kba=sqrt(nba/2)
bascaled<-t(scale(t(as.matrix(ba))))
colnames(bascaled)<-colnames(regcounts)
rownames(bascaled)<-rownames(ba)
bas<-cor(t(bascaled),method="pearson")
bap<-cor(t(bascaled),method="spearman")
bak<-cor(t(bascaled),method="kendall")


#    C L U S T E R I N G
#       Distance matrix from raw data
badist<-dist(ba)
#       Distance matrix from correlation matrix
bacordist<-dist(1-bacor)
#       Distance matrix from scaled data
bascaledist<-dist(bascaled)


# Hierarchical clustering
baclust<-hclust(badist,method="centroid")
basub<-cutree(baclust,k=kba)
#write.table(basub,file="circos/data/interactions/hierarchical_butyrate.raw",sep="\t",quote=FALSE,col.names=FALSE)
sils<-silhouette(basub,badist)
# K means clustering
km<-kmeans(badist,kba,iter.max=24,nstart=8)
write.table(km$cluster,file="circos/data/interactions/hierarchical_butyrate.raw",sep="\t",quote=FALSE,col.names=FALSE)
sils<-silhouette(km$cluster,badist)
# DBSCAN (make iteration?)
dr<-dbscan(badist,0.15,method="dist",showplot=TRUE)
#write.table(data.frame(rownames(bascaled)[dr$scaled > 0],dr$cluster[dr$cluster > 0]),file="circos/data/interactions/hierarchical_butyrate.raw",sep="\t",quote=FALSE,col.names=FALSE)

#       OPTICS - DONT FORGET TO RUN THE ALGORITHM


write.table(ba,file="clustering/raw/raw.csv",sep=",",quote=F,row.names=T,col.names=F)
write.table(bap,file="clustering/pearson/pearson.csv",sep=",",quote=F,row.names=T,col.names=F)
write.table(bas,file="clustering/spearman/spearman.csv",sep=",",quote=F,row.names=T,col.names=F)
write.table(bak,file="clustering/kendall/kendall.csv",sep=",",quote=F,row.names=T,col.names=F)
write.table(bascaled,file="clustering/scaled/scaled.csv",sep=",",quote=F,row.names=T,col.names=F)



clusters<-read.table(file="clustering/summary.csv",sep=",",col.names=c("subdir","data","minpts","reach","maxima","clust","dist","area","silhouette","davies","dunn","intracluster-dissimilarity","min-inter","max-inter"))
clusteres<-read.table(file="clustering/final/optics-clustering.csv",sep=",",col.names=c("id","cluster"))

x<-regcounts[as.character(clusteres[clusteres$cluster == "1",]$id),]
y<-length(x[,1])
x<-melt(x)
x$tex<-rep(c(rep("None",2*x),rep("TEX",x),rep("None",x),rep("TEX",x),rep("None",5*x)),3)
x$rep<-rep(c(rep("A",6*x),rep("B",4*x)),3)
x$time<-rep(c(rep("15",x),rep("150",x),rep("270",2*x),rep("75",2*x),rep("15",x),rep("150",x),rep("270",x),rep("75",x)),3)
x$stress<-c(rep("BA",10*x),rep("BuOH",10*x),rep("NS",10*x))
ggplot(x)+geom_smooth(aes(x=time,y=value,colour=stress))




# Partitioning around medoids
pam(badist,kba,diss=TRUE)


# Butanol
buoh<-regcounts[rownames(subset(stress[[1]],subset=stress[[1]]$padj < a & abs(stress[[1]]$log2FoldChange) > l)),]
nbuoh=length(rownames(buoh))
kbuoh=sqrt(nbuoh/2)
buohscaled<-t(scale(t(as.matrix(buoh))))
colnames(buohscaled)<-colnames(regcounts)
rownames(buohscaled)<-rownames(buoh)
buohcor<-cor(t(buohscaled),method="pearson")
buohdist<-as.dist(1-buohcor)
buohclust<-hclust(buohdist,method="centroid")
buohsub<-cutree(buohclust,k=kbuoh)
#write.table(buohsub,file="circos/data/interactions/hierarchical_butanol.raw",sep="\t",quote=FALSE,col.names=FALSE)
# Normality test
apply(regcounts,2,shapiro.test)



