# SUMMARY.R
# Matt Ralston
# Updated: 5/19/2014
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
#library("cairoDevice")
#Cairo()

mydir="/home/mrals/ETP/"
setwd(mydir)
source("summary/functions.r")
#                         A l i g n m e n t    S u m m a r y


# Alignment numbers and rRNA removal
# NOTE: rename files as NS-270-A or similar to facilitate spliting by factor
align<-read.table("summary-alignment.txt",header=TRUE,sep="\t")
align$totalaligned <- align$conc_once + align$conc_mult + align$disc_once + align$mate_once + align$mate_mult + align$unpair_once + align$unpair_mult
paired<-read.table("summary-paired.txt",header=TRUE,sep="\t")
unpaired<-read.table("summary-unpaired.txt",header=TRUE,sep="\t")
summary<-cbind(paired[,2]+unpaired[,2],paired[,2],unpaired[,2],paired[,3],unpaired[,3],align[,c(2,16)])
colnames(summary)<-c("Total","Paired","Unpaired","Paired ","Unpaired ","rRNA-free", "Aligned")
summary<-melt(summary)

x<-cbind(tapply(summary$value,summary$variable,mean),tapply(summary$value,summary$variable,sd),c("0","Trimmed","Trimmed","rRNA-free","rRNA-free","rRNA-free","0"))
x<-data.frame(factor(rownames(x),levels=rownames(x),ordered=TRUE),as.numeric(x[,1]),as.numeric(x[,2]),factor(x[,3]))
colnames(x)<-c("Class","Mean","Std.dev.","pair")
p1<-ggplot(x,aes(x=Class,y=Mean))+geom_bar(stat="identity",aes(fill=factor(pair)))+scale_fill_manual(name="Legend",breaks=c("Trimmed","rRNA-free"),values=c("grey34","darkgreen","blue4"))+theme(axis.title.x=element_blank())+scale_y_continuous(breaks=(0:4)*10**7)+geom_errorbar(aes(ymin=Mean-Std.dev.,ymax=Mean+Std.dev.))+annotate("text",label=c("55\u00B111%","45\u00B15%","19\u00B14%","15\u00B12%","35\u00B17%","35\u00B17%"),x=c(2:7),y=rep(2e6,6),colour=c(rep("white",5),"white"))
ggsave("summary/images/alignment_summary.png",p1,height=4,width=10)
# Mapping quality
mapq<-read.table("summary/mapq_summary.txt",header=TRUE)[,(1:5)]
colnames(mapq)<-c("NS270A","NS30A","NS30B","NS75A","NS75B")
mapq<-melt(as.matrix(mapq))
p1<-ggplot(mapq)+ geom_boxplot(aes(x=Var2,y=value))
ggsave("summary/images/mapq_boxplot.png",width=6,height=3,plot=p1)





#                         A s s e m b l y    S u m m a r y
assemblies<-read.table("summary/assemblies.size.txt")
p1<-ggplot(assemblies)+geom_boxplot(aes(x=V1,y=V1))+scale_y_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l')+xlab("All transcripts")+ylab("Transcript length, bp")+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())
p2<-ggplot(data.frame(V1=assemblies[assemblies$V1 < 5000,]))+geom_boxplot(aes(x=V1,y=V1))+xlab("Transcripts < 5kb (93%)")+ylab("Transcript length, bp")+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())
p3<-ggplot(data.frame(V1=assemblies[assemblies$V1 < 2000,]))+geom_boxplot(aes(x=V1,y=V1))+xlab("Transcripts < 2kb (80%)")+ylab("Transcript length, bp")+theme(axis.ticks.x=element_blank(),axis.text.x=element_blank())
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave.wide("summary/images/assembly_summary.png",p)


#                       G E N E   C O V E R A G E   S U M M A R Y

files<-list.files(paste(mydir,"summary/coverage",sep="/"),pattern="*.avcov",full.names=T)
my.list<-list()
my.frame<-data.frame()
col<-c()
for (i in 1:length(files)) {
    cov<-read.table(files[i])
    cov<-cbind(cov,rowMeans(cov[,(2:101)]))
    #covdev<-read.table("summary/NS30A+.sdcov")
    colnames(cov)<-c("Gene_id",(1:100), "Avg")
    mcov<-melt(cov[,2:101])
    colnames(mcov)<-c("percent","avg")
# Partitions data into quartiles
    mq1<-melt(cov[cov$Avg > summary(cov$Avg)[1] & cov$Avg < summary(cov$Avg)[2],(2:101)])
    mq2<-melt(cov[cov$Avg > summary(cov$Avg)[2] & cov$Avg < summary(cov$Avg)[3],(2:101)])
    mq3<-melt(cov[cov$Avg > summary(cov$Avg)[3] & cov$Avg < summary(cov$Avg)[5],(2:101)])
    mq4<-melt(cov[cov$Avg > summary(cov$Avg)[5] & cov$Avg < summary(cov$Avg)[6],(2:101)])
    
    
    q<-data.frame((1:100),t(apply(cov[,(2:101)],2,median.quartile)))
    colnames(q)<-c("percent","second","middle","fourth")
    q1<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[1] & cov$Avg < summary(cov$Avg)[2],(2:101)],2,median.quartile))))
    q2<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[2] & cov$Avg < summary(cov$Avg)[3],(2:101)],2,median.quartile))))
    q3<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[3] & cov$Avg < summary(cov$Avg)[5],(2:101)],2,median.quartile))))
    q4<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[5] & cov$Avg < summary(cov$Avg)[6],(2:101)],2,median.quartile))))
# Bar plot
    p1<-ggplot(q1,aes(x=V1,y=y))+geom_bar(stat="identity")+xlab("Percentage of gene")+ylab("Median Coverage of 1st quartile")+scale_x_discrete(breaks=pretty_breaks(n=10))
    p2<-ggplot(q2,aes(x=V1,y=y))+geom_bar(stat="identity")+xlab("Percentage of gene")+ylab("Median Coverage of 2nd quartile")+scale_x_discrete(breaks=pretty_breaks(n=10))
    p3<-ggplot(q3,aes(x=V1,y=y))+geom_bar(stat="identity")+xlab("Percentage of gene")+ylab("Median Coverage of 3rd quartile")+scale_x_discrete(breaks=pretty_breaks(n=10))
    p4<-ggplot(q4,aes(x=V1,y=y))+geom_bar(stat="identity")+xlab("Percentage of gene")+ylab("Median Coverage of 4th quartile")+scale_x_discrete(breaks=pretty_breaks(n=10))
    p<-arrangeGrob(p1,p2,p3,p4,ncol=2)
    ggsave.square(filename=paste(mydir,"summary/images/Quartile_cov_",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p)
    
    x<-as.matrix(cbind(as.numeric(mcov$variable), mcov$value))
    #hist(x[,2],n=200)
    
    
    # SCATTER
    p1<-ggplot(mcov,aes(percent,avg))+geom_point(alpha=0.2,size=1.5)+ylab("Coverage")+xlab("Percentage of gene")+scale_x_discrete(breaks=pretty_breaks(n=10))+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")+stat_summary(fun.y=median.quartile,geom="point",colour="darkred")
# Jitter
    p2<-ggplot(mcov,aes(x=percent,y=avg))+geom_jitter(alpha=0.2,size=1)+ylab("Coverage")+xlab("Percentage of gene")+scale_x_discrete(breaks=pretty_breaks(n=10))+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")+stat_summary(fun.y=median.quartile,geom="point",colour="red")
    
# Violin plot
    p3<-ggplot(mcov,aes(percent,avg))+geom_violin(fill="grey4")+stat_summary(fun.y=median.quartile,geom='point',colour="red")+ylab("Coverage")+xlab("Percentage of gene")+scale_x_discrete(breaks=pretty_breaks(n=10))+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")
    ggsave(paste("summary/images/scatter_cov",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p1)
    ggsave(paste("summary/images/jitter_cov",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p2)
    ggsave(paste("summary/images/violin_cov",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p3)
                                        # Half violin
#ggplot(data=mcov,aes(x=log10(avg)))+stat_density(aes(y=..density..))+scale_x_continuous(labels=math_format(10^.x))+facet_grid(. ~ percent)+coord_flip()
    my.list[[i]]<-mcov
    col[i]<-strsplit(tail(strsplit(files[i],"/")[[1]],n=1),"\\.")[[1]][1]
}

                     
for (i in 1:length(my.list)) {
    my.frame<-cbind.fill(my.frame,my.list[[i]]$avg)
}
colnames(my.frame)<-col
mframe<-melt(my.frame)
p1<-ggplot(mframe,aes(Var2,value))+geom_violin(fill="grey4")+stat_summary(fun.y=median.quartile,geom='point',colour="red")+ylab("Coverage")+xlab("Sample")+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")
ggsave("summary/images/Violin_coverage.png",p1,width=8,height=4)

#                        WHOLE   C o v e r a g e   S u m m a r y
# This shows how to generate a plot from coverage data
# One plot is created for each strand and plasmid


files<-list.files(paste(mydir,"Expression",sep="/"),pattern="*.cov",full.names=T)
for (i in 1:length(files)) {
  tmp<-read.table(files[i],header=FALSE)
  chrom<-tmp[tmp$V1 != 'gi|15004705|ref|NC_001988.2|',]
  plas<-tmp[tmp$V1 == 'gi|15004705|ref|NC_001988.2|',]
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


##############################################################################
#
#                                   D E S E Q 2
#
##############################################################################

directory<-"Expression"
sampleFiles<-list.files(paste(mydir,directory,sep="/"),pattern="*.3.counts",full.names=T)
sampleTime<-c(270,30,30,75,75)
sampleReplicate<-c(1,1,2,1,2)
#sampleTreatment<-c("NS","BuOH","Butyrate")
samples<-data.frame(sampleName=sampleFiles,fileName=sampleFiles,time=sampleTime,replicate=sampleReplicate)
samples<-samples[order(samples$time),]
row.names(samples)<-c(1,2,3,4,5)
htseq<-DESeqDataSetFromHTSeqCount(sampleTable=samples,directory='',design= ~ time + replicate)
#htseq<-DESeqDataSetFromHTSeqCount(sampleTable=samples,directory=directory,design= ~ time + treatment)
htseq$time<-factor(htseq$time,levels=c(30,75,270))
htseq$replicate<-factor(htseq$replicate,levels=c(1,2))
#htseq$treatment<-factor(htseq$treatment,levels=c("NS","BuOH","Butyrate"))
deseq<-DESeq(htseq)
rdeseq<-rlogTransformation(deseq,blind=TRUE)
#rdeseq<-rlog(deseq,blind=FALSE)


#                                  R A W   C O U N T S
cnt<-as.data.frame(counts(deseq,normalized=FALSE))
colnames(cnt)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
mcnt<-melt(cnt)
# Violin1
p1<-ggplot(mcnt)+geom_violin(aes(x=variable,y=value),trim=TRUE)+stat_summary(fun.y=median,fun.ymin=bottom.quartile,fun.ymax=top.quartile,geom='crossbar',aes(x=variable,y=value),col="red")+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l')
ggsave("summary/images/Counts_vertical.png",p1)
# Violin2
p2<-ggplot(mcnt)+geom_violin(aes(x=variable,y=value),trim=TRUE)+stat_summary(fun.y=median,fun.ymin=bottom.quartile,fun.ymax=top.quartile,geom='crossbar',aes(x=variable,y=value),col="red")+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+coord_flip()
ggsave("summary/images/Counts_horizontal.png",p2)
                                        # Histogram
p1<-ggplot(mcnt,aes(x=value, fill=variable))+geom_histogram(binwidth=.5,alpha=.2,position="identity")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," counts per gene")))+ylab("Frequency")+labs(fill="Legend")
ggsave("summary/images/Counts_histogram.png",p1)
# Histogram + Density estimate (probability
p2<-ggplot(mcnt,aes(x=value))+geom_histogram(aes(y=..density..,fill=variable),binwidth=.5,alpha=.1, position="identity")+geom_density(aes(colour=variable))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," counts per gene")))+ylab("Probability Density")+labs(fill="Legend")
ggsave("summary/images/Counts_density.png",p2)
#  C O R R E L A T I O N
p1<-ggplot(data=cnt,aes(x=NS30A +1,y=NS30B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," counts per gene: 30min, rep. A"))) + ylab(expression(paste(Log[10]," counts per gene: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS30A + 1) ~ log10(NS30B + 1), cnt)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS30A+1,y=NS30B+1,colour="best fit"),data=cnt,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Count Correlation: NS-30min.")
p2<-ggplot(data=cnt,aes(x=NS75A +1,y=NS75B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," counts per gene: 30min, rep. A"))) + ylab(expression(paste(Log[10]," counts per gene: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS75A + 1) ~ log10(NS75B + 1), cnt)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS75A+1,y=NS75B+1,colour="best fit"),data=cnt,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Count Correlation: NS-75min.")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave.wide("summary/images/Correlation_counts.png",p)



#                                 N O R M A L I Z E D     C O U N T S
normcounts<-as.data.frame(counts(deseq,normalized=TRUE))
colnames(normcounts)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
mnorm<-melt(normcounts)
# Histogram
p1<-ggplot(mnorm,aes(x=value, fill=variable))+geom_histogram(binwidth=.5,alpha=.2,position="identity")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," normalized counts per gene")))+ylab("Frequency")
ggsave("summary/images/Norm_histogram.png",p1)
# Histogram + Density estimate (probability)
p2<-ggplot(mnorm,aes(x=value))+geom_histogram(aes(y=..density..,fill=variable),binwidth=.5,alpha=.1, position="identity")+geom_density(aes(colour=variable))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," normalized counts per gene")))+ylab("Probability Density")
ggsave("summary/images/Norm_density.png",p2)
# Correlation of size-normalized counts
p1<-ggplot(data=normcounts,aes(x=NS30A +1,y=NS30B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS30A + 1) ~ log10(NS30B + 1), normcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS30A+1,y=NS30B+1,colour="best fit"),data=normcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Normalized Count Correlation: NS-30min.")
p2<-ggplot(data=normcounts,aes(x=NS75A +1,y=NS75B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS75A + 1) ~ log10(NS75B + 1), normcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS75A+1,y=NS75B+1,colour="best fit"),data=normcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Normalized Count Correlation: NS-75min.")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave.wide("summary/images/Correlation_norm_counts.png",p)


#                                 R E G U L A R I Z E D   C O U N T S
regcounts<-as.data.frame(assay(rdeseq))
colnames(regcounts)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
mreg<-melt(regcounts/log2(10))
# Histogram
p1<-ggplot(mreg,aes(x=(value+1)/log2(10), fill=variable))+geom_histogram(binwidth=.5,alpha=.1,position="identity")+scale_x_continuous(breaks=(0:5),labels=math_format(10^.x))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," regularized counts per gene")))+ylab("Frequency")
ggsave("summary/images/Regularized_histogram.png",p1)
# Histogram + Density estimate (probability)
p2<-ggplot(mreg,aes(x=(value+1)/log2(10)))+geom_histogram(aes(y=..density..,fill=variable),binwidth=.5,alpha=.1, position="identity")+geom_density(aes(colour=variable))+scale_x_continuous(breaks=(0:5),labels=math_format(10^.x))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," regularized counts per gene")))+ylab("Probability Density")
ggsave("summary/images/Regularized_density.png",p2)
# Correlation of regularized counts
p1<-ggplot(data=regcounts,aes(x=(NS30A+1)/log2(10),y=(NS30B+1)/log2(10))) + geom_point(alpha=0.5,size=2)+scale_y_continuous(labels=math_format(10^.x))+scale_x_continuous(labels=math_format(10^.x))+annotation_logticks(base=10,sides='bl') + xlab(expression(paste(Log[10]," regularized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," regularized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 1, y = 3, label = lm_eqn(lm(NS30A+1 ~ NS30B+1, regcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=(NS30A+1)/log2(10),y=(NS30B+1)/log2(10),colour="best fit"),data=regcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Regularized Count Correlation: NS-30min.")
p2<-ggplot(data=regcounts,aes(x=(NS75A+1)/log2(10),y=(NS75B+1)/log2(10))) + geom_point(alpha=0.5,size=2)+scale_y_continuous(labels=math_format(10^.x))+scale_x_continuous(labels=math_format(10^.x))+annotation_logticks(base=10,sides='bl') + xlab(expression(paste(Log[10]," regularized counts: 75min, rep. A"))) + ylab(expression(paste(Log[10]," regularized counts: 75min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 1, y = 3, label = lm_eqn(lm(NS75A+1 ~ NS75B+1, regcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=(NS75A+1)/log2(10),y=(NS75B+1)/log2(10),colour="best fit"),data=regcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Regularized Count Correlation: NS-75min.")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave.wide("summary/images/Correlation_reg_counts.png",p)



#                                 R E S U L T S
res<-results(deseq)
res30v270<-results(deseq,contrast=c("time","270","30"))
res30v75<-results(deseq,contrast=c("time","75","30"))
res75v270<-results(deseq,contrast=c("time","270","75"))


#                                   M A    p l o t s

p1<-ggplot(as.data.frame(res30v75),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(low="red",high="blue")+ylab(expression(paste(log[2],"75min/30min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 75min/30min")
p2<-ggplot(as.data.frame(res75v270),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(low="red",high="blue")+ylab(expression(paste(log[2],"270min/75min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 270min/75min")
p3<-ggplot(as.data.frame(res30v270),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(low="red",high="blue")+ylab(expression(paste(log[2],"270min/30min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 270min/30min")
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave.wide("summary/images/MA-plots.png",p)


#                             V A R I A N C E   +    D I S P E R S I O N    E S T I M A T I O N  + F I T T I N G
countmodel<-lm(NS30A~NS30B,data=log10(cnt+1))
normmodel<-lm(NS30A~NS30B,data=log10(normcounts+1))
regmodel<-lm(NS30A~NS30B,data=(regcounts/log2(10)))
len<-length(residuals(countmodel))

countmodel<-data.frame(rep("Counts",len),fitted(countmodel),residuals(countmodel))
normmodel<-data.frame(rep("Norm. counts",len),fitted(normmodel),residuals(normmodel))
regmodel<-data.frame(rep("Reg. counts",len),fitted(regmodel),residuals(regmodel))
colnames(countmodel)<-c("V1","V2","V3")
colnames(normmodel)<-c("V1","V2","V3")
colnames(regmodel)<-c("V1","V2","V3")
resid<-rbind(countmodel,normmodel,regmodel)
colnames(resid)<-c("Variable","Fitted","Residuals")
#                 R E S I D U A L    P L O T
# Fit becomes more homoscedastic
p1<-ggplot(resid)+geom_point(aes(x=Fitted,y=Residuals,colour=Variable),alpha=0.7,size=2)+scale_colour_manual("Legend",values=c("Black","#56B4E9","red"),labels=c("Dispersion",expression(paste("Max.",italic(' a posteriori'))),"Model fit"))
ggsave('summary/images/Residual_plot.png',p1)

#                 D I S P E R S I O N   P L O T
mseq<-as.data.frame(mcols(deseq))
disp<-melt(as.data.frame(mcols(deseq))[,c(1,4,9,5)],id="baseMean")
p1<-ggplot(disp)+geom_point(aes(x=log10(baseMean),y=log10(value),colour=variable),alpha=0.7,size=2) + scale_colour_manual("Legend",values=c("black","#56B4E9","red"),labels=c("Dispersion",expression(paste("Max.",italic(' a posteriori'))),"Model fit")) + xlab(expression(paste(Log[10]," regularized expression"))) + ylab(expression(paste(Log[10]," dispersion")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+theme(legend.justification=c(0,0),legend.position=c(0,0))+guides(colour=guide_legend(override.aes=list(size=5)))
ggsave("summary/images/Gene_dispersion.png",plot=p1,width=5,height=3)


#                             V O L C A N O    P L O T

p1<-ggplot(as.data.frame(res30v75),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"75min/30min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 75 min vs 30 min")
p2<-ggplot(as.data.frame(res75v270),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"270min/75min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 270 min vs 75 min")
p3<-ggplot(as.data.frame(res30v270),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"270min/30min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 270 min vs 30 min")
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave.wide("summary/images/Volcano_plot.png",p)

#      P r i n c i p a l      C o m p o n e n t      A n a l y s i s

png("summary/images/PCA.png")
plotPCA(rdeseq,intgroup=c("time"))
dev.off()

# Heatmap
# select genes
select <- order(rowMeans(counts(deseq,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# Standard read counts
heatmap.2(counts(deseq,normalized=TRUE)[select,], col = hmcol,Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10,6))
# Variance stabilizing transformation (similar to log read counts)
vsd<-assay(varianceStabilizingTransformation(deseq))[,c(2,3,4,5,1)]
colnames(vsd)<-c("30A","30B","75A","75B","270")
heatmap.2(vsd[select,], col = hmcol, Rowv = FALSE, Colv = FALSE, scale="none", dendrogram="none", trace="none", margin=c(10, 6))
# Differentially expressed genes (in this case res came from the comparison NS30 vs NS270)
tmp<-res[complete.cases(res),]
degenes<-counts(deseq,normalized=TRUE)[rownames(tmp[tmp$padj < 0.05,]),]
# Heatmap with dendrogram of DE genes and optional subset
heatmap.2(vsd[rownames(tmp[tmp$padj < 0.05,]),], col = hmcol, Colv = FALSE, scale="none", dendrogram="row", trace="none", key=FALSE, margin=c(10, 6))
#heatmap.2(vsd[rownames(tmp[tmp$padj < 0.05,])[1:50],], col = hmcol, Colv = FALSE, scale="none", dendrogram="row", trace="none", key=FALSE, margin=c(10, 6))


