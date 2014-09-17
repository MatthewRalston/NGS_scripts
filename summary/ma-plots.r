# MA-plots.R
# Copyright 2014 Matt Ralston
# Updated 7/26/2014
#
# This R script contains code to generate MA-plots from DESeq2
# Differential expression analyses

library('gtools')

times<-c("15","75","150","270")
stresses<-c("BA","BuOH","NS")
texes<-c("Normal","TEX")

time1=list()
time2=list()

buohtime=list()
batime=list()

# Time only
for (i in 2:4) {
    #  (ti vs t1) 
    t<-results(deseq,contrast=list(paste('time',times[i],sep=''),paste('time',times[1],sep='')))
    p1<-ggplot(as.data.frame(t),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(bquote(Log[2] ~ .(paste("(",times[i],"min/15min)",sep=''))))+xlab(expression(paste(Log[10],"regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-3:3),labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col='blue')+geom_hline(yintercept=c(-0.5,0.5),col='dodgerblue',lwd=2)+ggtitle(paste("MA-plot: ",times[i]," min / 15 min", sep=''))
    ggsave(paste("summary/images/maplots//time/T",i,"vT1.png",sep=''),p1)
    time1[[i]]=t
    t<-as.data.frame(t)
    l<-length(rownames(t))
    cndtn<-rep(paste("T",i),l)
    cntrl<-rep("T1",l)
    t<-data.frame(ID=rownames(t),Logexp=log10(t$baseMean),foldchange=t$log2FoldChange,pval=t$padj,condition=cndtn,control=cntrl)
    write.table(t,file="summary/maplot.csv",quote=F,sep=",",row.names=F,col.names=F,append=T)
}
for (i in 3:4) {
    # ti vs ti-1
    t<-results(deseq,contrast=list(paste('time',times[i],sep=''),paste('time',times[i-1],sep='')))
    p1<-ggplot(as.data.frame(t),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(bquote(Log[2] ~ .(paste("(",times[i],"min/",times[i-1],"min)",sep=''))))+xlab(expression(paste(Log[10],"regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-3:3),labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col='blue')+geom_hline(yintercept=c(-0.5,0.5),col='dodgerblue',lwd=2)+ggtitle(paste("MA-plot: ",times[i]," min / ",times[i-1]," min", sep=''))
    ggsave(paste("summary/images/maplots/time/T",i,'vT',i-1,".png",sep=''),p1)
    time2[[i]]=t
    t<-as.data.frame(t)
    l<-length(rownames(t))
    cndtn<-rep(paste("T",i),l)
    cntrl<-rep(paste("T",i-1),l)
    t<-data.frame(ID=rownames(t),Logexp=log10(t$baseMean),foldchange=t$log2FoldChange,pval=t$padj,condition=cndtn,control=cntrl)
    write.table(t,file="summary/maplot.csv",quote=F,sep=",",row.names=F,col.names=F,append=T)
}

# Stress only (BA vs NS) (BuOH vs NS) (BA vs BuOH)
ba<-results(deseq,contrast=list('condBA','condNS'))
buoh<-results(deseq,contrast=list('condBuOH','condNS'))
p1<-ggplot(as.data.frame(buoh),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(expression(paste(Log[2],"(BuOH/NS)")))+xlab(expression(paste(Log[10],"regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-3:3),labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col='blue')+geom_hline(yintercept=c(-0.5,0.5),col='dodgerblue',lwd=2)+ggtitle("MA-plot: BuOH / NS")
p2<-ggplot(as.data.frame(ba),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(expression(paste(Log[10],"(p)")),low="red",high="blue")+ylab(expression(paste(Log[2],"(BA/NS)")))+xlab(expression(paste(Log[10],"regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-3:3),labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col='blue')+geom_hline(yintercept=c(-0.5,0.5),col='dodgerblue',lwd=2)+ggtitle("MA-plot: BA / NS")
ggsave("summary/images/maplots/stress/BuOHvsNS.png",p1)
ggsave("summary/images/maplots/stress/BAvsNS.png",p2)
stress=list(buoh=buoh,butyrate=ba)

ba<-as.data.frame(ba)
buoh<-as.data.frame(buoh)
l<-length(rownames(ba))
butanol<-rep("Butanol",l)
butyrate<-rep("Butyrate",l)
cntrl<-rep("NS",l)
buoh<-data.frame(ID=rownames(buoh),Logexp=log10(buoh$baseMean),foldchange=buoh$log2FoldChange,pval=buoh$padj,condition=butanol,control=cntrl)
ba<-data.frame(ID=rownames(ba),Logexp=log10(ba$baseMean),foldchange=ba$log2FoldChange,pval=ba$padj,condition=butyrate,control=cntrl)
write.table(ba,file="summary/maplot.csv",quote=F,sep=",",row.names=F,col.names=F,append=T)
write.table(buoh,file="summary/maplot.csv",quote=F,sep=",",row.names=F,col.names=F,append=T)

# Stress over time (BA vs NS) (BuOH vs NS) at each time
for (i in 1:length(times)) {
    butyrate<-results(deseq,contrast=list(paste('time',times[i],'.cond',stresses[1],sep=''),paste('time',times[i],'.cond',stresses[3],sep='')))
    butanol<-results(deseq,contrast=list(paste('time',times[i],'.cond',stresses[2],sep=''),paste('time',times[i],'.cond',stresses[3],sep='')))
    x="(padj)"
    if (length(levels(factor(butanol$padj))) == 1) {
        butanol$padj<-butanol$pvalue
        x="(p)"
    }
    p1<-ggplot(as.data.frame(butanol),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(bquote(Log[10] ~ .(eval(x))),low="red",high="blue")+ylab(bquote(Log[2] ~ .(paste("(BuOH",times[i],"min/NS",times[i],"min)"))))+xlab(expression(paste(Log[10],"regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-3:3),labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col='blue')+geom_hline(yintercept=c(-0.5,0.5),col='dodgerblue',lwd=2)+ggtitle(paste("MA-plot: BuOH ",times[i],"min / NS ",times[i],"min",sep=''))
    p2<-ggplot(as.data.frame(butyrate),aes(log10(baseMean),log2FoldChange))+geom_point(aes(colour=log10(padj)),alpha=0.7)+scale_colour_gradient(bquote(Log[10] ~ .('(padj)')),low="red",high="blue")+ylab(bquote(Log[2] ~ .(paste("(BA",times[i],"min/NS",times[i],"min)"))))+xlab(expression(paste(Log[10],"regularized counts")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=(-3:3),labels=math_format(2^.x))+annotation_logticks(sides='b')+geom_hline(yintercept=0,col='blue')+geom_hline(yintercept=c(-0.5,0.5),col='dodgerblue',lwd=2)+ggtitle(paste("MA-plot: BA ",times[i],"min / NS ",times[i],"min",sep=''))
    ggsave(paste("summary/images/maplots/time_stress/BuOHvNS_",times[i],".png",sep=''),p1)
    ggsave(paste("summary/images/maplots/time_stress/BAvNS_",times[i],".png",sep=''),p2)
    batime[[i]]=butyrate
    buohtime[[i]]=butanol


    ba<-as.data.frame(butyrate)
    buoh<-as.data.frame(butanol)
    l<-length(rownames(buoh))
    butanol<-rep(paste("Butanol-",times[i]),l)
    butyrate<-rep(paste("Butyrate-",times[i]),l)
    cntrl<-rep(paste("NS",times[i]),l)
    buoh<-data.frame(ID=rownames(buoh),Logexp=log10(buoh$baseMean),foldchange=buoh$log2FoldChange,pval=buoh$padj,condition=butanol,control=cntrl)
    ba<-data.frame(ID=rownames(ba),Logexp=log10(ba$baseMean),foldchange=ba$log2FoldChange,pval=ba$padj,condition=butyrate,control=cntrl)
    write.table(ba,file="summary/maplot.csv",quote=F,sep=",",row.names=F,col.names=F,append=T)
    write.table(buoh,file="summary/maplot.csv",quote=F,sep=",",row.names=F,col.names=F,append=T)
}
