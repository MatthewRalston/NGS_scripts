library('ggplot2')
library('scales')
library("edgeR")
library("DESeq2")
library("RColorBrewer")
library("gplots")
require("gridExtra")
library("reshape2")
library("cairoDevice")
Cairo()

mydir="/home/mrals/ETP/"
setwd(mydir)

#                         A l i g n m e n t    S u m m a r y
# Mapping quality

mapq<-read.table("summary/mapq_summary.txt",header=TRUE)[,(1:5)]
colnames(mapq)<-c("NS270A","NS30A","NS30B","NS75A","NS75B")
mapq<-melt(as.matrix(mapq))
p1<-ggplot(mapq)+ geom_boxplot(aes(x=Var2,y=value))


# Alignment numbers and rRNA removal
# NOTE: rename files as NS-270-A or similar to facilitate spliting by factor
align<-read.table("summary-alignment.txt",header=TRUE,sep="\t")
align$totalaligned <- align$conc_once + align$conc_mult + align$disc_once + align$mate_once + align$mate_mult + align$unpair_once + align$unpair_mult
paired<-read.table("summary-paired.txt",header=TRUE,sep="\t")
unpaired<-read.table("summary-unpaired.txt",header=TRUE,sep="\t")
summary<-cbind(paired[,2]+unpaired[,2],paired[,2],unpaired[,2],paired[,3],unpaired[,3],align[,c(2,16)])
colnames(summary)<-c("Total reads","Paired reads","Unpaired reads","rRNA-free paired reads","rRNA-free unpaired reads","Total rRNA-free reads", "Total aligned reads")
summary<-cbind(melt(summary), rep(c(0,1,1,2,2,0,0),each=5))
colnames(summary)[3]<-"pair"
ggplot(summary) + geom_boxplot(aes(x=variable,y=value,fill=factor(pair)))+scale_fill_manual(name="Legend",values=c("white","grey","blue"))


#                         A s s e m b l y    S u m m a r y
assemblies<-read.table("summary/assemblies.size.txt")
boxplot(assemblies, ylab="Transcript size, bp")
boxplot(assemblies[assemblies$V1 < 5000,],ylab="Transcript size, bp")



#                           C o v e r a g e   S u m m a r y
# This shows how to generate a plot from coverage data
# One plot is created for each strand and plasmid

system('rm coverage/*.png')
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

#                         C o u n t    C o r r e l a t i o n s
count<-read.table("countsummary.txt",header=TRUE)
count$id <- c(as.factor(count$Gene_id))
ggplot(data=count,aes(NS30A,NS30B))+geom_point(alpha=0.3)+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='bl')+xlab(expression(paste(Log[10]," counts : 30 minutes, rep. A")))+ylab(expression(paste(Log[10]," counts : 30 minutes, rep. B")))
abline(0,1,col="red")

ggplot(data=count,aes(NS75A,NS75B))+geom_point(alpha=0.3)+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='bl')+xlab(expression(paste(Log[10]," counts : 30 minutes, rep. A")))+ylab(expression(paste(Log[10]," counts : 30 minutes, rep. B")))
abline(0,1,col="red")
NS30<-count[(count$NS30A!=0 & count$NS30B!=0),]
NS75<-count[(count$NS75A!=0 & count$NS75B!=0),]
summary(NS30.lm<-lm(log(NS30A)~log(NS30B),data=NS30))
summary(NS75.lm<-lm(log(NS75A)~log(NS75B),data=NS75))

ggplot(count[,c(2,3,4,5,6)],aes(x=c("30A","30B","75A","75B","270"),y=c(NS270A NS30A NS30B NS75A NS75B)))+geom_boxplot()


# Geometric mean normalization (Cufflinks) 
# Check the vignette for the quality control of this technique
# (size factors, dispersion estimates, PCA analysis)
geom<-read.table("Cuffquant/genes.count_table.geometric", header=TRUE)
# obnoxiously highly expressed genes
geom[(geom$NS30_0 > 500 | geom$NS30_1 > 500 | geom$NS75_0 > 500 | geom$NS75_1 > 500 | geom$NS270_0 > 500 ),]
ggplot(data=geom,aes(x=NS30_0,y=NS30_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS30_0) ~ log10(NS30_1), geom[geom$NS30_0 > 0.00 & geom$NS30_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS30_0),y=log10(NS30_1),colour="best fit"),data=geom[geom$NS30_0 > 0.00 & geom$NS30_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
ggplot(data=geom,aes(x=NS75_0,y=NS75_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10]," Cufflinks normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS75_0) ~ log10(NS75_1), geom[geom$NS75_0 > 0.00 & geom$NS75_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS75_0),y=log10(NS75_1),colour="best fit"),data=geom[geom$NS75_0 > 0.00 & geom$NS75_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
summary(NS30.geom<-lm(log(NS30_0)~log(NS30_1),data=geom))
summary(NS75.geom<-lm(log(NS75_0)~log(NS75_1),data=geom))


#FPKM Normalized expression
fpkm<-read.table("Cuffquant/genes.count_table.fpkm",header=TRUE)
#*****
png("Correlation_fpkm_counts.png",width=1440,height=400)
p1<-ggplot(data=fpkm,aes(x=NS30_0,y=NS30_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS30_0) ~ log10(NS30_1), fpkm[fpkm$NS30_0 > 0.00 & fpkm$NS30_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS30_0),y=log10(NS30_1),colour="best fit"),data=fpkm[fpkm$NS30_0 > 0.00 & fpkm$NS30_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p2<-ggplot(data=fpkm,aes(x=NS75_0,y=NS75_1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. A"))) + ylab(expression(paste(Log[10],"(FPKM) normalized expression: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 5, y = 1000, label = lm_eqn(lm(log10(NS75_0) ~ log10(NS75_1), fpkm[fpkm$NS75_0 > 0.00 & fpkm$NS75_1 > 0.00,])), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=log10(NS75_0),y=log10(NS75_1),colour="best fit"),data=fpkm[fpkm$NS75_0 > 0.00 & fpkm$NS75_1 > 0.00,],method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
grid.arrange(p1,p2,ncol=2)
dev.off()
summary(NS30.fpkm<-lm(log(NS30_0)~log(NS30_1),data=fpkm))
summary(NS75.fpkm<-lm(log(NS75_0)~log(NS75_1),data=fpkm))



#         D E S E Q,  N o r m a l i z a t i o n,  R e g u l a r i z a t i o n
directory<-"counts"
files<-list.files(paste(mydir,"Expression",sep="/"),pattern="*.3.counts",full.names=T)
sampleTime<-c(270,30,30,75,75)
#sampleTreatment<-c("NS","BuOH","Butyrate")
samples<-data.frame(sampleName=sampleFiles,fileName=sampleFiles,time=sampleTime)
samples<-samples[order(samples$time),]
row.names(samples)<-c(1,2,3,4,5)
htseq<-DESeqDataSetFromHTSeqCount(sampleTable=samples,directory=directory,design= ~ time)
#htseq<-DESeqDataSetFromHTSeqCount(sampleTable=samples,directory=directory,design= ~ time + treatment)
htseq$time<-factor(htseq$time,levels=c(30,75,270))
#htseq$treatment<-factor(htseq$treatment,levels=c("NS","BuOH","Butyrate"))
deseq<-DESeq(htseq)
rdeseq<-rlog(deseq,blind=TRUE)
#rdeseq<-rlog(deseq,blind=FALSE)
#                                 C O U N T S
normcounts<-as.data.frame(counts(deseq,normalized=TRUE))
colnames(normcounts)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
cnt<-as.data.frame(counts(deseq,normalized=FALSE))
colnames(cnt)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
regcounts<-as.data.frame(assay(rdeseq))
colnames(regcounts)<-c("NS30A","NS30B","NS75A","NS75B","NS270A")
#                                 R E S U L T S
res<-results(deseq)
res30v270<-results(deseq,contrast=c("time","270","30"))
res30v75<-results(deseq,contrast=c("time","75","30"))
res75v270<-results(deseq,contrast=c("time","270","75"))
#                             C O R R E L A T I O N
# plot of counts
png("Correlation_counts.png",width=1440,height=400)
p1<-ggplot(data=cnt,aes(x=NS30A +1,y=NS30B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," counts per gene: 30min, rep. A"))) + ylab(expression(paste(Log[10]," counts per gene: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS30A + 1) ~ log10(NS30B + 1), cnt)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS30A+1,y=NS30B+1,colour="best fit"),data=cnt,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p2<-ggplot(data=cnt,aes(x=NS75A +1,y=NS75B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," counts per gene: 30min, rep. A"))) + ylab(expression(paste(Log[10]," counts per gene: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS75A + 1) ~ log10(NS75B + 1), cnt)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS75A+1,y=NS75B+1,colour="best fit"),data=cnt,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
grid.arrange(p1,p2,ncol=2)
dev.off()
# Plot of size-normalized counts
png("Correlation_norm_counts.png",width=1440,height=400)
p1<-ggplot(data=normcounts,aes(x=NS30A +1,y=NS30B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS30A + 1) ~ log10(NS30B + 1), normcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS30A+1,y=NS30B+1,colour="best fit"),data=normcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p2<-ggplot(data=normcounts,aes(x=NS75A +1,y=NS75B + 1)) + geom_point(alpha=0.5,size=2)+scale_y_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+scale_x_log10(breaks=10**(0:4),labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides='l') + xlab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," geom. normalized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 20, y = 2000, label = lm_eqn(lm(log10(NS75A + 1) ~ log10(NS75B + 1), normcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=NS75A+1,y=NS75B+1,colour="best fit"),data=normcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
grid.arrange(p1,p2,ncol=2)
dev.off()
# Plot of regularized counts
png("Correlation_reg_counts.png",width=1440,height=400)
p1<-ggplot(data=regcounts,aes(x=(NS30A+1)/log2(10),y=(NS30B+1)/log2(10))) + geom_point(alpha=0.5,size=2)+scale_y_continuous(labels=math_format(10^.x))+scale_x_continuous(labels=math_format(10^.x))+annotation_logticks(base=10,sides='bl') + xlab(expression(paste(Log[10]," regularized counts: 30min, rep. A"))) + ylab(expression(paste(Log[10]," regularized counts: 30min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 1, y = 3, label = lm_eqn(lm(NS30A+1 ~ NS30B+1, regcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=(NS30A+1)/log2(10),y=(NS30B+1)/log2(10),colour="best fit"),data=regcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
p2<-ggplot(data=regcounts,aes(x=(NS75A+1)/log2(10),y=(NS75B+1)/log2(10))) + geom_point(alpha=0.5,size=2)+scale_y_continuous(labels=math_format(10^.x))+scale_x_continuous(labels=math_format(10^.x))+annotation_logticks(base=10,sides='bl') + xlab(expression(paste(Log[10]," regularized counts: 75min, rep. A"))) + ylab(expression(paste(Log[10]," regularized counts: 75min, rep. B")))+geom_abline(intercept=0,slope=1,size=1,col="blue")+annotate("text", x = 1, y = 3, label = lm_eqn(lm(NS75A+1 ~ NS75B+1, regcounts)), colour="black", size = 5, parse=TRUE)+stat_smooth(aes(x=(NS75A+1)/log2(10),y=(NS75B+1)/log2(10),colour="best fit"),data=regcounts,method="lm",size=1,na.rm=TRUE,fullrange=TRUE)+scale_color_discrete("")
grid.arrange(p1,p2,ncol=2)
dev.off()
#                                   M A    p l o t s
png("MA-plots.png",width=1440,height=400)
p1<-ggplot(as.data.frame(res30v75),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(low="red",high="blue")+ylab(expression(paste(log[2],"75min/30min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 75min/30min")
p2<-ggplot(as.data.frame(res75v270),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(low="red",high="blue")+ylab(expression(paste(log[2],"270min/75min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 270min/75min")
p3<-ggplot(as.data.frame(res30v270),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(low="red",high="blue")+ylab(expression(paste(log[2],"270min/30min")))+xlab(expression(paste(Log[10]," regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col="blue") + geom_hline(yintercept=seq(-0.5,0.5,1),col="dodgerblue",lwd=2)+ggtitle("MAplot: 270min/30min")
grid.arrange(p1,p2,p3,ncol=3)
dev.off()


#                                 D I S P E R S I O N
disp<-melt(as.data.frame(mcols(deseq))[,c(1,4,9,5)],id="baseMean")
p1<-ggplot(disp)+geom_point(aes(x=log10(baseMean),y=log10(value),colour=variable),alpha=0.7,size=2) + scale_colour_manual("Legend",values=c("black","#56B4E9","red"),labels=c("Dispersion",expression(paste("Max.",italic(' a posteriori'))),"Model fit")) + xlab(expression(paste(Log[10]," regularized expression"))) + ylab(expression(paste(Log[10]," dispersion")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+annotation_logticks(sides='b')+theme(legend.justification=c(0,0),legend.position=c(0,0))+guides(colour=guide_legend(override.aes=list(size=5)))
ggsave(filename="Gene_dispersion.png",plot=p1,width=5,height=3)


#                             V O L C A N O    P L O T
png("Volcano_plot.png",width=1440,height=400)
p1<-ggplot(as.data.frame(res30v75),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"75min/30min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 75 min vs 30 min")
p2<-ggplot(as.data.frame(res75v270),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"270min/75min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 270 min vs 75 min")
p3<-ggplot(as.data.frame(res30v270),aes(log2FoldChange,-log10(padj)))+geom_point(aes(colour=-log10(padj),aplha=0.7))+scale_colour_gradient(low="red",high="green")+xlab(expression(paste(log[2],"270min/30min")))+ylab(expression(-log[10]("p.value")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(labels=math_format(10^.x))+geom_vline(xintercept=0,col="blue")+geom_vline(xintercept=seq(-0.5,0.5,1),col="dodgerblue")+ggtitle("Volcano plot: 270 min vs 30 min")
grid.arrange(p1,p2,p3,ncol=3)
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



#                                   F U N C T I O N S

# Function- Multiplot

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



lm_eqn = function(m) {
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}


