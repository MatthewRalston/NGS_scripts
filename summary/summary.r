# SUMMARY.R
# Copyright 2014 Matt Ralston
# Updated: 7/23/2014
#
# This R script contains code to generate summary and diagnostic plots
# for RNA-seq data. Included are correlations, dispersion estimates,
# heatmaps, per gene coverage, and much more.
library('ggplot2')
library('scales')
library("edgeR")
library("DESeq2")
library("RColorBrewer")
library("gplots")
require("gridExtra")
library("reshape2")
library("plotly")
#require('rgl')
#library("MASS")
#library("cairoDevice")

options(max.print=1000)
#Cairo()
plot.ly <- plotly(username="MatthewRalston",key="yd9yj8vmx4")
mydir="/home/mrals/Final/"
setwd(mydir)
source("summary/functions.r")
tss_window=100
colnum=10
rownum=10
#                         A l i g n m e n t    S u m m a r y

source("summary/alignmentstats.r")

#            D i f f e r e n t i a l    E x p r e s s i o n
source("summary/diff_exp.r")

# Mapping quality
mapq<-read.table("summary/mapq_summary.txt",header=TRUE)[,(1:5)]
colnames(mapq)<-c("NS270A","NS30A","NS30B","NS75A","NS75B")
mapq<-melt(as.matrix(mapq))
p1<-ggplot(mapq)+ geom_boxplot(aes(x=Var2,y=value))
ggsave("summary/images/mapq_boxplot.png",width=6,height=3,plot=p1)





#                         A s s e m b l y    S u m m a r y
assemblies<-read.table("summary/assemblies.size.txt")
p1<-ggplot(assemblies)+geom_histogram(aes(x=V1))+scale_x_log10(breaks=10**(-1:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='b')+ylab("Counts")+xlab("Transcript size, bp")
p2<-ggplot(data.frame(V1=assemblies[assemblies$V1 < 5000,]))+geom_boxplot(aes(x=V1,y=V1))+xlab("Transcripts < 5kb (93%)")+ylab("Transcript size, bp")+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())
p3<-ggplot(data.frame(V1=assemblies[assemblies$V1 < 2000,]))+geom_boxplot(aes(x=V1,y=V1))+xlab("Transcripts < 2kb (80%)")+ylab("Transcript size, bp")+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave.wide("summary/images/assembly_summary.png",p)

p<-ggplot(assemblies)+geom_violin(aes(x="test",y=V1))+scale_y_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l')
ggsave("summary/images/assembly_violin.png",p)


ggsave("summary/images/Violin_coverage.png",p1,width=8,height=4)


#                        WHOLE   C o v e r a g e   S u m m a r y
# This shows how to generate a plot from coverage data
# One plot is created for each strand and plasmid


files<-list.files(paste(mydir,"Expression",sep="/"),pattern="*.cov",full.names=T)
for (i in 1:length(files)) {
  tmp<-read.table(files[i],header=FALSE)
  chrom<-tmp[tmp$V1 != 'NC_001988.2',]
  plas<-tmp[tmp$V1 == 'NC_001988.2',]
  png(filename=paste(strsplit(files[i],split="\\.")[[1]][1],"chrom.png",sep=""),width=1440,height=400)
  
  p1<-ggplot(chrom,aes(x=V2,y=V3)) + geom_point(alpha=0.08,size=1)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + geom_rug(alpha=0.08,sides='r') + xlab("Coordinates, bp") + ylab("Coverage")
  p2<-ggplot(chrom,aes(x='log10 counts',y=V3),axes=FALSE) + geom_boxplot()+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l')
  multiplot(p1,p2+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()),layout=matrix(c(1,1,1,1,2),nrow=1))
  dev.off()
  
  
  png(filename=paste(strsplit(files[i],split="\\.")[[1]][1],"plas.png",sep=""),width=1440,height=400)
  p1<-ggplot(plas,aes(x=V2,y=V3)) + geom_point(alpha=0.08,size=1)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + geom_rug(alpha=0.08,sides='r') + xlab("Coordinates, bp") + ylab("Coverage")
  p2<-ggplot(plas,aes(x='log10 counts',y=V3),axes=FALSE) + geom_boxplot()+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l')
  multiplot(p1,p2+theme(axis.text.x=element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank()),layout=matrix(c(1,1,1,1,2),nrow=1))
  dev.off()
  
}




# Geometric mean normalization (Cufflinks) 
# Check the vignette for the quality control of this technique
# (size factors, dispersion estimates, PCA analysis)
geom<-read.table("Cuffquant/genes.count_table.geometric", header=TRUE)
# obnoxiously highly expressed genes
geom[(geom$NS30_0 > 500 | geom$NS30_1 > 500 | geom$NS75_0 > 500 | geom$NS75_1 > 500 | geom$NS270_0 > 500 ),]
p1<-ggplot(data=geom,aes(x=NS30_0,y=NS30_1)) + geom_point(alpha=0.7,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS30_0) ~ log10(NS30_1), geom[geom$NS30_0 > 0.00 & geom$NS30_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS30_0),y=log10(NS30_1),colour="best fit"),data=geom[geom$NS30_0 > 0.00 & geom$NS30_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p2<-ggplot(data=geom,aes(x=NS75_0,y=NS75_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS75_0) ~ log10(NS75_1), geom[geom$NS75_0 > 0.00 & geom$NS75_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS75_0),y=log10(NS75_1),colour="best fit"),data=geom[geom$NS75_0 > 0.00 & geom$NS75_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave.wide("summary/images/Cufflinks_correlation.png",p)


#summary(NS30.geom<-lm(log(NS30_0)~log(NS30_1),data=geom))
#summary(NS75.geom<-lm(log(NS75_0)~log(NS75_1),data=geom))


#FPKM Normalized expression
fpkm<-read.table("Cuffquant/genes.count_table.fpkm",header=TRUE)
#*****

p1<-ggplot(data=fpkm,aes(x=NS30_0,y=NS30_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS30_0) ~ log10(NS30_1), fpkm[fpkm$NS30_0 > 0.00 & fpkm$NS30_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS30_0),y=log10(NS30_1),colour="best fit"),data=fpkm[fpkm$NS30_0 > 0.00 & fpkm$NS30_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p2<-ggplot(data=fpkm,aes(x=NS75_0,y=NS75_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS75_0) ~ log10(NS75_1), fpkm[fpkm$NS75_0 > 0.00 & fpkm$NS75_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS75_0),y=log10(NS75_1),colour="best fit"),data=fpkm[fpkm$NS75_0 > 0.00 & fpkm$NS75_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave.wide("summary/images/FPKM_correlation.png",p)

#summary(NS30.fpkm<-lm(log(NS30_0)~log(NS30_1),data=fpkm))
#summary(NS75.fpkm<-lm(log(NS75_0)~log(NS75_1),data=fpkm))


# Store sets of upregulated and downregulated genes according to the contrast
c.frame<-data.frame()
cutoff<-0.05
for (i in 1:length(contlist)) {
    up<-(contlist[[i]]$padj < cutoff & contlist[[i]]$log2FoldChange > 0)
    up<-replace(up, is.na(up),FALSE)
    down<-(contlist[[i]]$padj < cutoff & contlist[[i]]$log2FoldChange < 0)
    down<-replace(down, is.na(down),FALSE)
    x<-replace(replace(rep("NONE",length(res$padj)),up,"UP"),down,"DOWN")
    c.frame<-cbind.fill(c.frame,x)
}
colnames(c.frame)<-contrasts
rownames(c.frame)<-rownames(res)
uplist<-list()
downlist<-list()
for (i in 1:length(contlist)) {
    uplist[[i]]<-rownames(c.frame[c.frame[,i] == "UP",])
    downlist[[i]]<-rownames(c.frame[c.frame[,i] == "DOWN",])
}
names(uplist)<-contrasts
names(downlist)<-contrasts
# Create a Venn diagram of the set intersections
venn.diagram(uplist,main="Upregulated genes",'summary/images/Venn_up.tiff')
venn.diagram(downlist,main="Downregulated genes",'summary/images/Venn_down.tiff')

# Create list of set intersections #http://stat.ethz.ch/R-manual/R-devel/library/base/html/sets.html
uptime<-intersect(uplist$`30v270`,uplist$`30v75`)
downtime<-intersect(downlist$`30v270`,downlist$`30v75`)
# INTERACTIVE ONLY: Print a list of upregulated genes for analysis in DAVID
#cat(uplist$`30v75`,sep="\n")




##############################################################################
#
#         T R A N S C R I P T I O N     S T A R T    S I T E S
#
###############################################################################
# This code is used to create ggplots of coverage near the tss


files<-list.files(paste(mydir,"summary/coverage",sep="/"),pattern="*.tss",full.names=T)

for (i in 1:length(files)) {
    tsses <- read.table(files[i],header=T)
    colnames(tsses) <- c("Gene_id",(1:tss_window-1))
    # 3D plots
    mtss <- melt(tsses[,2:tss_window])
    mtss.kde <- kde2d(as.numeric(mtss[,1]),log10(x[,2]+1),n=500)
    col1<-heat.colors(length(mtss.kde$z))[rank(mtss.kde$z)]
    #axes3d(c('x','y','z'),col='black',labels=FALSE)
    #title3d('','','x','y','z')
    #persp3d(x=mtss.kde,col=col1)
    persp(mtss.kde,phi=15,theta=0,shade=.3,border=NA,xlab="Position (bp)",ylab="Log10(Coverage)",zlab="Density")
    # Jitter plot of coverage versus position
    ggplot(mtss,aes(x=variable,y=value+1))+geom_jitter(alpha=0.3,size=1)+xlab("Position (bp)")+ylab("Coverage")+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")+stat_summary(fun.y=median.quartile,geom="point",colour="red")
    # Histogram of coverage across the first base
    
    pltList=list()
    for (i in 1:10){
        for (j in 1:10){
            p<-ggplot(tsses,aes(x=tsses[,(i-1)*10+j]))+geom_histogram()+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10],"(Coverage)")))+ggtitle(paste("Coverage at base ", (i-1)*10+j))
           plotList=c(plotList,list(p))
        }
    }
    do.call(grid.arrange, list(plotList,colnum))
}
