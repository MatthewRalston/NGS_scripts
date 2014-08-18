# SUMMARY.R
# Copyright 2014 Matt Ralston
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
library("plotly")
#require('rgl')
#library("MASS")
#library("cairoDevice")ssh

#Cairo()
plot.ly <- plotly(username="MatthewRalston",key="yd9yj8vmx4")
mydir="/home/mrals/Final/"
setwd(mydir)
source("summary/functions.r")
tss_window=100
colnum=10
rownum=10
#                         A l i g n m e n t    S u m m a r y


# Alignment numbers and rRNA removal
# NOTE: rename files as NS-270-A or similar to facilitate spliting by factor
talign<-read.table("summary-alignment.txt",header=TRUE,sep="\t")
align<-talign[!grepl("TEX",talign$ID),]
texalign<-talign[grepl("TEX",talign$ID),]
align$totalaligned <- align$conc_once + align$conc_mult + align$disc_once + align$mate_once + align$mate_mult + align$unpair_once + align$unpair_mult
texalign$totalaligned <- texalign$conc_once + texalign$conc_mult + texalign$disc_once + texalign$mate_once + texalign$mate_mult + texalign$unpair_once + texalign$unpair_mult
tpaired<-read.table("summary-paired.txt",header=TRUE,sep="\t")
paired<-tpaired[!grepl("TEX",tpaired$ID),]
texpaired<-tpaired[grepl("TEX",tpaired$ID),]
tunpaired<-read.table("summary-unpaired.txt",header=TRUE,sep="\t")
unpaired<-tunpaired[!grepl("TEX",tunpaired$ID),]
texunpaired<-tunpaired[grepl("TEX",tunpaired$ID),]

summary<-cbind(paired[,2]+unpaired[,2],paired[,2],unpaired[,2],paired[,3],unpaired[,3],align[,c(2,16)],rep("Normal",length(paired[,2])))
colnames(summary)<-c("Total","Trimmed Paired","Trimmed Unpaired","Paired","Unpaired","rRNA-free", "Aligned","Tex")

texsummary<-cbind(texpaired[,2]+texunpaired[,2],texpaired[,2],texunpaired[,2],texpaired[,3],texunpaired[,3],texalign[,c(2,16)],rep("TEX",length(texpaired[,2])))
colnames(texsummary)<-c("Total","Trimmed Paired","Trimmed Unpaired","Paired","Unpaired","rRNA-free", "Aligned","Tex")
texsummary<-melt(texsummary)
summary<-melt(summary)


testsum<-data.frame(factor(c(rep("Raw reads",24),rep("Trimmed",48),rep("rRNA-free",72),rep("Total Aligned",24)),ordered=FALSE),summary$variable,c(rep("Total",24),rep("Paired",24),rep("Unpaired",24),rep("Paired",24),rep("Unpaired",24),rep("Total",48)),summary$Tex,summary$value)

texsum<-data.frame(factor(c(rep("Raw reads",6),rep("Trimmed",12),rep("rRNA-free",18),rep("Total Aligned",6)),ordered=FALSE),texsummary$variable,c(rep("Total",6),rep("Paired",6),rep("Unpaired",6),rep("Paired",6),rep("Unpaired",6),rep("Total",12)),texsummary$Tex,texsummary$value)
x<-c("one","two","three","four","five")
colnames(testsum)<-x
colnames(texsum)<-x
test<-data.frame(rbind(testsum,texsum))
colnames(test)<-c("grid","group","mycolor","tex","value")
test$grid<-factor(test$grid,levels(test$grid)[c(1,4,2,3)])
test$mycolor<-factor(test$mycolor,c("Paired","Unpaired","Total"))

text<-data.frame(grid=c(rep("Trimmed",4),rep("rRNA-free",6),rep("Total Aligned",2)),mycolor=c(rep(c("#003300","#0000FF"),3),"#999999","#003300","#0000FF",rep("#999999",3)), tex=c(rep("Normal",2),rep("TEX",2),rep("Normal",3),rep("TEX",3),"Normal","TEX"), value=c(4e6,0,8.5e6,1e6,1.43e6,0.9e6,6.3e6,3e6,0,10e6,1.4e6,4.5e6), lab=c("90.1%", "9.8%","84.6%","15.3%", "34.3%", "3.5%","37.8%","41%","6.8%","47.9%","39.7%", "50.2%"))

p1<-ggplot(text,aes(x=tex,y=value))+geom_violin(aes(fill=mycolor,scale="width"),scale="width")+facet_grid(~grid)+theme(axis.title.x=element_blank())+scale_fill_manual(name="Legend",values=c("darkgreen","blue4","grey34"))+geom_text(data=text,aes(x=tex,y=value,label=lab,colour=mycolor,size=0.7))+scale_colour_manual(values=c("blue4","darkgreen","grey34"),guide="none")+scale_size(guide="none")+ylab("Reads per Library")+scale_y_continuous(breaks=seq(0,12,2)*10**6)



ggsave("summary/images/alignment_summary.png",p1,height=4,width=10)
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


#                       G E N E   C O V E R A G E   S U M M A R Y

#** Update summary.r to calculate/plot:
#*** standard deviation of coverage per gene
#*** coefficient of variation of coverage per gene (sigma/mu)
#*** Plot avg coverage vs gene length
#*** Plot coefficient of variation vs gene length
#*** Plot avg coverage vs GC content
#*** Plot coefficient of variation vs GC

bias<-data.frame(a=numeric(),b=numeric(),c=numeric(),d=numeric())
files<-list.files(paste(mydir,"summary/coverage",sep="/"),pattern="*.avcov",full.names=T)
my.list<-list()
my.frame<-data.frame()
col<-c()
for (i in 1:length(files)) {
    cov<-read.table(files[i],header=T)
    cov$Avg<-rowMeans(cov[,(4:103)])
    cov$sd<-apply(cov[,(4:103)],1,sd)
    cov$cov<-cov$Avg*cov$sd
    colnames(cov)<-c("Gene_id","length","gc", (1:100), "Avg","sd","cov")
    mcov<-melt(cov[,4:103])
    colnames(mcov)<-c("percent","avg")
    # Check for correlation between Avg coverage, coefficient of variation x length, gc content
    x<-log10(cov$Avg + 1)
    y<-log10(cov$cov + 1)
    z<-c(summary(lm(x~cov$length))$adj.r.squared,summary(lm(x~cov$gc))$adj.r.squared, summary(lm(y~cov$length))$adj.r.squared, summary(lm(y~cov$gc))$adj.r.squared)
    bias<-rbind(bias,z)
    
    
# Partitions data into quartiles
    mq1<-melt(cov[cov$Avg > summary(cov$Avg)[1] & cov$Avg < summary(cov$Avg)[2],(4:103)])
    mq2<-melt(cov[cov$Avg > summary(cov$Avg)[2] & cov$Avg < summary(cov$Avg)[3],(4:103)])
    mq3<-melt(cov[cov$Avg > summary(cov$Avg)[3] & cov$Avg < summary(cov$Avg)[5],(4:103)])
    mq4<-melt(cov[cov$Avg > summary(cov$Avg)[5] & cov$Avg < summary(cov$Avg)[6],(4:103)])
    
    
    q<-data.frame((1:100),t(apply(cov[,(4:103)],2,median.quartile)))
    colnames(q)<-c("percent","second","middle","fourth")
    q1<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[1] & cov$Avg < summary(cov$Avg)[2],(4:103)],2,median.quartile))))
    q2<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[2] & cov$Avg < summary(cov$Avg)[3],(4:103)],2,median.quartile))))
    q3<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[3] & cov$Avg < summary(cov$Avg)[5],(4:103)],2,median.quartile))))
    q4<-as.data.frame(cbind((1:100),t(apply(cov[cov$Avg > summary(cov$Avg)[5] & cov$Avg < summary(cov$Avg)[6],(4:103)],2,median.quartile))))
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
colnames(bias)<-c("avg.vs.len","av.vs.gc","cov.vs.len","cov.vs.gc")
              
for (i in 1:length(my.list)) {
    my.frame<-cbind.fill(my.frame,my.list[[i]]$avg)
}
colnames(my.frame)<-col
mframe<-melt(my.frame)
p1<-ggplot(mframe,aes(Var2,value))+geom_violin(fill="grey20")+stat_summary(fun.y=median.quartile,geom='point',colour="red")+ylab("Coverage")+xlab("Sample")+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")
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


##############################################################################
#
#                                   D E S E Q 2
#
##############################################################################

directory<-"Expression/counts"
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
cnt<-as.data.frame(counts(deseq,normalized=FALSE))
colnames(cnt)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
normcounts<-as.data.frame(counts(deseq,normalized=TRUE))
colnames(normcounts)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
regcounts<-as.data.frame(assay(rdeseq))
colnames(regcounts)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")


#                                  R A W   C O U N T S
mcnt<-melt(cnt)
# Violin1
p1<-ggplot(mcnt,aes(x=variable,y=value))+geom_violin(fill='grey20',trim=TRUE)+stat_summary(fun.y=median.quartile,geom='point',col="red")+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l')+theme(axis.title.x=element_blank())+ylab("Counts per gene")
ggsave("summary/images/Counts_vertical.png",p1)
# Violin2
p2<-ggplot(mcnt,aes(x=variable,y=value))+geom_violin(fill='grey20',trim=TRUE)+stat_summary(fun.y=median,fun.ymin=bottom.quartile,fun.ymax=top.quartile,geom='crossbar',aes(x=variable,y=value),col="red")+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+xlab("Sample")+ylab("Counts per gene")+theme(axis.title.x=element_blank())+coord_flip()
ggsave("summary/images/Counts_horizontal.png",p2)
# Jitter
p3<-ggplot(mcnt,aes(x=variable,y=value))+geom_jitter()+scale_y_log10(breaks=10**(0:4),labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='l')+stat_summary(fun.y=median,fun.ymax=top.quartile,fun.ymin=bottom.quartile,geom='crossbar',colour='red')+ylab("Counts")+theme(axis.title.x=element_blank())
ggsave("summary/images/Counts_jitter.png",p3)

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
mnorm<-melt(normcounts)
# Histogram
p1<-ggplot(mnorm,aes(x=value, fill=variable))+geom_histogram(binwidth=.5,alpha=.2,position="identity")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," normalized counts per gene")))+ylab("Frequency")
ggsave("summary/images/Norm_histogram.png",p1)
# Histogram + Density estimate (probability)
p2<-ggplot(mnorm,aes(x=value))+geom_histogram(aes(y=..density..,fill=variable),binwidth=.5,alpha=.1, position="identity")+geom_density(aes(colour=variable))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," normalized counts per gene")))+ylab("Probability Density")
ggsave("summary/images/Norm_density.png",p2)
# Jitter 
p3<-ggplot(mnorm,aes(x=variable,y=value))+geom_jitter()+scale_y_log10(breaks=10**(0:4),labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='l')+stat_summary(fun.y=median,fun.ymax=top.quartile,fun.ymin=bottom.quartile,geom='crossbar',colour='red')+ylab("Normalized Counts")+theme(axis.title.x=element_blank())
ggsave("summary/images/Norm_jitter.png",p3)
# Correlation of size-normalized counts
p1<-ggplot(data=normcounts,aes(x=NS30A +1,y=NS30B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS30A + 1) ~ log10(NS30B + 1), normcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS30A+1,y=NS30B+1,colour="best fit"),data=normcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Normalized Count Correlation: NS-30min.")
p2<-ggplot(data=normcounts,aes(x=NS75A +1,y=NS75B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS75A + 1) ~ log10(NS75B + 1), normcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS75A+1,y=NS75B+1,colour="best fit"),data=normcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")+ggtitle("Normalized Count Correlation: NS-75min.")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave.wide("summary/images/Correlation_norm_counts.png",p)


#                                 R E G U L A R I Z E D   C O U N T S
mreg<-melt((regcounts+1)/log2(10))
# Histogram
p1<-ggplot(mreg,aes(x=value, fill=variable))+geom_histogram(binwidth=.5,alpha=.1,position="identity")+scale_x_continuous(breaks=(0:5),labels=math_format(10^.x))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," regularized counts per gene")))+ylab("Frequency")
ggsave("summary/images/Regularized_histogram.png",p1)
# Histogram + Density estimate (probability)
p2<-ggplot(mreg,aes(x=value))+geom_histogram(aes(y=..density..,fill=variable),binwidth=.5,alpha=.1, position="identity")+geom_density(aes(colour=variable))+scale_x_continuous(breaks=(0:5),labels=math_format(10^.x))+annotation_logticks(base=10,sides="b")+xlab(expression(paste(Log[10]," regularized counts per gene")))+ylab("Probability Density")
ggsave("summary/images/Regularized_density.png",p2)
# Jitter
p3<-ggplot(mreg,aes(x=variable,y=value))+geom_jitter()+scale_y_continuous(breaks=(0:4),labels=math_format(10^.x))+annotation_logticks(base=10,sides='l')+stat_summary(fun.y=median,fun.ymax=top.quartile,fun.ymin=bottom.quartile,geom='crossbar',colour='red')+ylab("Regularized Counts")+theme(axis.title.x=element_blank())
ggsave("summary/images/regularized_jitter.png",p3)
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
contlist<-list(res30v270,res75v270,res30v75)
contrasts<-c("30v270","75v270","30v75")

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

#                                   M A    p l o t s

p1<-ggplot(as.data.frame(res30v75),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(expression(paste(log[2],"75min/30min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 75min/30min")
ggsave("summary/images/MA_30_75.png",p1)
p2<-ggplot(as.data.frame(res75v270),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(expression(paste(log[2],"270min/75min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 270min/75min")
p3<-ggplot(as.data.frame(res30v270),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(expression(paste(log[2],"270min/30min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 270min/30min")
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave("summary/images/MA_30_270.png",p3)
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
p1<-ggplot(resid)+geom_point(aes(x=Fitted,y=Residuals,colour=Variable),alpha=0.7,size=2)+scale_colour_manual("Legend",values=c("Black","#56B4E9","red"),labels=c("Raw counts","Normalized counts",expression(paste("Max.",italic(' a posteriori')))))+xlab("Gene Expression")+scale_x_continuous(breaks=-1:4,labels=math_format(10^.x))+ylab("Residual")
ggsave('summary/images/Residual_plot.png',p1)

#                 D I S P E R S I O N   P L O T
mseq<-as.data.frame(mcols(deseq))
# Variance vs mean
p1<-ggplot(mseq,aes(x=baseMean,y=sqrt(baseVar)/baseMean))+geom_point()+scale_x_log10()+annotation_logticks(base=10,sides='b')+ylab("Coeff. of Var.")+xlab("Mean")
ggsave("summary/images/COV_mean.png",p1)
disp<-melt(as.data.frame(mcols(deseq))[,c(1,4,9,5)],id="baseMean")
# Dispersion plot
p1<-ggplot(disp)+geom_point(aes(x=log10(baseMean),y=log10(value),colour=variable),alpha=0.7,size=1.5) + scale_colour_manual("Legend",values=c("black","#56B4E9","red"),labels=c("Dispersion",expression(paste("Max.",italic(' a posteriori'))),"Model fit")) + xlab(expression(paste(Log[10]," regularized expression"))) + ylab(expression(paste(Log[10]," dispersion")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-7:2),labels=math_format(10^.x))+annotation_logticks(sides='b')+theme(legend.justification=c(0,0),legend.position=c(0,0))+guides(colour=guide_legend(override.aes=list(size=5)))

ggsave("summary/images/Gene_dispersion.png",plot=p1,width=5,height=3)
p1
response <- plotly$plotly(data


#                             V O L C A N O    P L O T

p1<-ggplot(as.data.frame(res30v75),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"75min/30min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 75 min vs 30 min")
p2<-ggplot(as.data.frame(res75v270),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"270min/75min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 270 min vs 75 min")
p3<-ggplot(as.data.frame(res30v270),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"270min/30min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 270 min vs 30 min")
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave.wide("summary/images/Volcano_plot.png",p)

#      P r i n c i p a l      C o m p o n e n t      A n a l y s i s

cormat<-cor(regcounts)
covmat<-cov(regcounts)
p1<-ggplot(melt(cormat))+geom_tile(aes(Var1,Var2,fill=value))
p2<-ggplot(melt(covmat))+geom_tile(aes(Var1,Var2,fill=value))
p<-arrangeGrob(p1,p2,ncol=2)
ggsave("summary/images/Correl_covar_heatmap.png",plot=p,width=4,height=2)
scree<-data.frame(variable=c("PC1","PC2","PC3","PC4","PC5"),value=prcomp(cormat)$sdev**2)
p1<-ggplot(scree)+geom_bar(aes(variable,value),stat='identity')+ylab("Variance")+theme(axis.title.x=element_blank())
ggsave("summary/images/pca_scree.png",p1)
v<-c("NS30A","NS30B","NS75A","NS75B","NS270")
mycols<-c("#990000","#33FF66","#009933","#3300FF","#6600CC")
corpca<-prcomp(regcounts)
melted<-melt(corpca$rotation)
p1<-ggplot(melted)+geom_bar(aes(x=Var1,y=value,fill=Var1),stat='identity')+facet_wrap(~Var2)
ggsave.square("summary/images/pca_barplot.png",p1)

scores<-data.frame(v,prcomp(cormat)$rotation)
p1<-ggplot(scores)+geom_point(aes(x=PC1,y=PC2,colour=v),stat='identity',size=5)+scale_colour_manual("Legend",values=mycols)+geom_segment(aes(x=0,y=0,xend=PC1*0.75,yend=PC2*0.75,colour=v),arrow=arrow(length=unit(0.3,"cm")))+labs("Legend")+xlab("PC:1")+ylab("PC:2")
p2<-ggplot(scores)+geom_point(aes(x=PC1,y=PC3,colour=v),stat='identity',size=5)+scale_colour_manual("Legend",values=mycols)+geom_segment(aes(x=0,y=0,xend=PC1*0.75,yend=PC3*0.75,colour=v),arrow=arrow(length=unit(0.3,"cm")))+labs("Legend")+xlab("PC:1")+ylab("PC:3")
p3<-ggplot(scores)+geom_point(aes(x=PC2,y=PC3,colour=v),stat='identity',size=5)+scale_colour_manual("Legend",values=mycols)+geom_segment(aes(x=0,y=0,xend=PC2*0.75,yend=PC3*0.75,colour=v),arrow=arrow(length=unit(0.3,"cm")))+labs("Legend")+xlab("PC:2")+ylab("PC:3")
ggsave("summary/images/PCA.png",p1)
p<-arrangeGrob(p1,p2,p3,ncol=3)
ggsave.wide("summary/images/pca_three.png",p)


require('rgl')
y<-c(0,0,0)
x<-scores[,2:4]
x<-rbind(y,x[1,],y,x[2,],y,x[3,],y,x[4,],y,x[5,])
rgl.open()
rgl.bg(color="white")
z<-x[1:4,]
rgl.lines(z[,1],z[,2],z[,3],color="blue")
z<-x[5:8,]
rgl.lines(z[,1],z[,2],z[,3],color="green")
z<-x[9:10,]
rgl.lines(z[,1],z[,2],z[,3],color="red")
axes3d(c('x','y','z'),col='black',labels=FALSE)
title3d('','','x','y','z')
n<-c(paste(rep("00",9),1:9,sep=''),paste(rep("0",89),10:99,sep=''),as.character(100:360))
for (i in 1:360) {
    rgl.viewpoint(i,0, interactive=F)
    rgl.snapshot(paste("summary/pca/file",n[i],".png",sep=""))
}
system('convert -delay 1 summary/pca/*.png summary/images/pca.gif')
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
