# AlignmentStats.R
# Copyright 2014 Matt Ralston
# Updated: 7/26/14
#
# This R script contains code to generate summary plots for the alignment
# of Fastq reads to the C_ac genome.
library('ggplot2')
library('scales')
library('cairoDevice')
library('gridExtra')
options(max.print=1000)
options(scipen=999)


mydir="/home/mrals/Final"
setwd(mydir)


# Coverage vector plot
names<-read.table("Terry-excercise/colnames.txt",sep=",")
m<-seq(1,59,2)
p<-seq(2,60,2)
minus<-names[m,1]
plus<-names[p,1]

cov<-read.table("Terry-excercise/CA_P0151-0176.cov",sep=",")
p1<-





#    A l i g n m e n t   S u m m a r y   V i o l i n   P l o t


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


sum<-data.frame(factor(c(rep("Raw reads",24),rep("Trimmed",48),rep("rRNA-free",72),rep("Total Aligned",24)),ordered=FALSE),summary$variable,c(rep("Total",24),rep("Paired",24),rep("Unpaired",24),rep("Paired",24),rep("Unpaired",24),rep("Total",48)),summary$Tex,summary$value)

texsum<-data.frame(factor(c(rep("Raw reads",6),rep("Trimmed",12),rep("rRNA-free",18),rep("Total Aligned",6)),ordered=FALSE),texsummary$variable,c(rep("Total",6),rep("Paired",6),rep("Unpaired",6),rep("Paired",6),rep("Unpaired",6),rep("Total",12)),texsummary$Tex,texsummary$value)
x<-c("one","two","three","four","five")
colnames(sum)<-x
colnames(texsum)<-x
final<-data.frame(rbind(sum,texsum))
colnames(final)<-c("grid","group","mycolor","tex","value")
final$grid<-factor(final$grid,levels(final$grid)[c(1,3,4,2)])
final$mycolor<-factor(final$mycolor,c("Paired","Unpaired","Total"))


# Percentages (total of normal, tex, are 100%
# Trimmed (Ratio of avg # of paired reads to avg total # of reads
# mean(paired$total)/(mean(paired$total)+mean(unpaired$total))
# rRNA-free
# total reads
# mean(align$total)/(mean(paired$total)+mean(unpaired$total))
# paired reads
# mean(align$paired)/(mean(paired$total)+mean(unpaired$total))
mycolors<-c(rep(c("#0000FF","#003300"),3),"#999999","#0000FF","#003300",rep("#999999",3))
values<-c(14e6,10e6,14e6,10e6,4e6,3e6,5e6,4e6,3e6,5e6,5e6,5e6)
labels<-c("89.5%","10.5","88.6%","11.4%","35.1%","3.8%","38.9%","43.4%","6.9%","50.3%","40.7%", "52.7%")
text<-data.frame(grid=c(rep("Trimmed",4),rep("rRNA-free",6),rep("Total Aligned",2)),mycolor=mycolors, tex=c(rep("Normal",2),rep("TEX",2),rep("Normal",3),rep("TEX",3),"Normal","TEX"), value=values, lab=labels)

p1<-ggplot(final,aes(x=tex,y=value))+geom_violin(aes(fill=mycolor,scale="width"),scale="width")+facet_grid(~grid)+theme(axis.title.x=element_blank())+scale_fill_manual(name="Legend",values=c("darkgreen","blue4","grey34"))+geom_text(data=text,aes(x=tex,y=value,label=lab,colour=mycolor,size=0.7),position=position_dodge(0.9))+scale_colour_manual(values=c("darkgreen","blue4","grey34"),guide="none")+scale_size(guide="none")+ylab("Reads per Library")+scale_y_continuous(breaks=seq(0,30,2)*10**6)+stat_summary(fun.y=median.quartile,fun.ymin=bottom.quartile,fun.ymax=top.quartile,geom='point',col='red',position=position_dodge(0.9),aes(fill=mycolor))



ggsave("summary/images/alignment_summary.png",p1,height=4,width=10)


#                      C o v e r a g e   S u m m a r y



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
    ggsave.square(filename=paste(mydir,"summary/images/coverage/Quartile_cov_",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p)
    
    x<-as.matrix(cbind(as.numeric(mcov$variable), mcov$value))
    #hist(x[,2],n=200)
    
    
    # SCATTER
    p1<-ggplot(mcov,aes(percent,avg))+geom_point(alpha=0.2,size=1.5)+ylab("Coverage")+xlab("Percentage of gene")+scale_x_discrete(breaks=pretty_breaks(n=10))+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")+stat_summary(fun.y=median.quartile,geom="point",colour="darkred")
# Jitter
    p2<-ggplot(mcov,aes(x=percent,y=avg))+geom_jitter(alpha=0.2,size=1)+ylab("Coverage")+xlab("Percentage of gene")+scale_x_discrete(breaks=pretty_breaks(n=10))+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")+stat_summary(fun.y=median.quartile,geom="point",colour="red")
    
# Violin plot
    p3<-ggplot(mcov,aes(percent,avg))+geom_violin(fill="grey4")+stat_summary(fun.y=median.quartile,geom='point',colour="red")+ylab("Coverage")+xlab("Percentage of gene")+scale_x_discrete(breaks=pretty_breaks(n=10))+scale_y_log10(breaks=10**(-1:4),limits=c(0.1,1e04),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="l")
    ggsave(paste("summary/images/coverage/scatter_cov",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p1)
    ggsave(paste("summary/images/coverage/jitter_cov",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p2)
    ggsave(paste("summary/images/coverage/violin_cov",tail(strsplit(files[i],"/")[[1]],n=1),".png",sep=""),p3)
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
ggsave("summary/images/coverage/summary_violin.png",p1)
