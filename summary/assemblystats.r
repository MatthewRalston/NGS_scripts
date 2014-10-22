# ASSEMBLYSTATS.R
# Copyright 2014 Matt Ralston
# Updated: 10/5/14
#
# This is the R script to produce distributions of transcript length,etc for the assemblies.
library('ggplot2')
library('scales')
library('cairoDevice')
library('gridExtra')
options(max.print=1000)
options(scipen=999)


mydir="/home/mrals/Final"
setwd(mydir)


# TOTAL
utrfive<-read.table("Trinity_ref/stats/five.len",header=F)
utrthree<-read.table("Trinity_ref/stats/three.len",header=F)
igr<-read.table("Trinity_ref/stats/igr.len",header=F)
refigr<-read.table("Trinity_ref/stats/refigr.len",header=F)
opsize<-read.table("Trinity_ref/stats/operonsize.txt",header=F)
oplen<-read.table("Trinity_ref/stats/operon.len",header=F)
novel<-read.table("Trinity_ref/stats/novel.len",header=F)
standard<-read.table("Trinity_ref/stats/standard.len",header=F)
refcds<-read.table("Trinity_ref/stats/refcds.len",header=F)
intra<-read.table("Trinity_ref/stats/intrautr.len",header=F)
# PAIRED
putrfive<-read.table("Trinity_paired/stats/five.len",header=F)
putrthree<-read.table("Trinity_paired/stats/three.len",header=F)
pigr<-read.table("Trinity_paired/stats/igr.len",header=F)
prefigr<-read.table("Trinity_paired/stats/refigr.len",header=F)
popsize<-read.table("Trinity_paired/stats/operonsize.txt",header=F)
poplen<-read.table("Trinity_paired/stats/operon.len",header=F)
pnovel<-read.table("Trinity_paired/stats/novel.len",header=F)
pstandard<-read.table("Trinity_paired/stats/standard.len",header=F)
prefcds<-read.table("Trinity_paired/stats/refcds.len",header=F)
pintra<-read.table("Trinity_paired/stats/intrautr.len",header=F)
ptotal<-(read.table("Trinity_paired/stats/total.len",header=F))$V1
# With annotation
futrfive<-read.table("annotation/stats/five.len",header=F)
futrthree<-read.table("annotation/stats/three.len",header=F)
figr<-read.table("annotation/stats/igr.len",header=F)
frefigr<-read.table("annotation/stats/refigr.len",header=F)
fopsize<-read.table("annotation/stats/operonsize.txt",header=F)
foplen<-read.table("annotation/stats/operon.len",header=F)
fnovel<-read.table("annotation/stats/novel.len",header=F)
fstandard<-read.table("annotation/stats/standard.len",header=F)
frefcds<-read.table("annotation/stats/refcds.len",header=F)
fintra<-read.table("annotation/stats/intrautr.len",header=F)
ftotal<-c(fnovel$V2,fstandard$V2)

paredesoplength<-read.table("summary/paredes-operon.len",header=F)
paredesopsize<-read.table("summary/paredes-operon-size.txt",header=F)
paredesintrautr<-read.table("summary/paredes-intra-utr.len",header=F)
paredesigr<-read.table("summary/paredes-igrs.len",header=F)
total<-c(novel$V1,standard$V1)




# Total data
utr<-rbind(data.frame(variable=rep("5`",length(utrfive$V1)),value=utrfive$V1),data.frame(variable=rep("3`",length(utrthree$V1)),value=utrthree$V1))

intra<-rbind(data.frame(variable=rep("Ref. Intra-UTR",length(paredesintrautr$V1)),value=paredesintrautr$V1),data.frame(variable=rep("Intra-UTR",length(intra$V1)),value=intra$V1))

igrs<-rbind(data.frame(variable=rep("Trinity",length(igr$V1)),value=igr$V1),data.frame(variable=rep("Reference",length(refigr$V1)),value=refigr$V1))
lens<-rbind(data.frame(variable=rep("Novel",length(novel$V1)),value=novel$V1),data.frame(variable=rep("Standard",length(standard$V1)),value=standard$V1),data.frame(variable=rep("Total",length(total)),value=total))

totallens<-rbind(data.frame(variable=rep("Novel",length(novel$V1)),value=novel$V1),data.frame(variable=rep("Standard",length(standard$V1)),value=standard$V1),data.frame(variable=rep("Operon",length(oplen$V1)),value=oplen$V1),data.frame(variable=rep("CDS",length(refcds$V1)),value=refcds$V1),data.frame(variable=rep("Ref. Operon",length(paredesoplength$V1)),value=paredesoplength$V1),data.frame(variable=rep("Total",length(total)),value=total))

oplength<-rbind(data.frame(variable=rep("CDS",length(refcds$V1)),value=refcds$V1),data.frame(variable=rep("Ref. Operon",length(paredesoplength$V1)),value=paredesoplength$V1),data.frame(variable=rep("Standard",length(standard$V1)),value=standard$V1),data.frame(variable=rep("Operon",length(oplen$V1)),value=oplen$V1))

operonsize<-rbind(data.frame(variable=rep("Standard",length(opsize$V1)),value=opsize$V1),data.frame(variable=rep("Ref. Operons",length(paredesopsize$V1[paredesopsize$V1 > 1])),value=paredesopsize$V1[paredesopsize$V1 > 1]))
# Total Graphing
p1<-ggplot(utr,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("UTR length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p2<-ggplot(intra,aes(x=value,fill=variable))+geom_density(alpha=0.3)+xlab("UTR length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')

p3<-ggplot(igrs,aes(x=value,fill=variable))+geom_density(alpha=0.3)+xlab("IGR length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')

p4<-ggplot(lens,aes(x=value,fill=variable))+geom_density(alpha=0.4)+xlab("Transcript length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p5<-ggplot(totallens,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("Feature length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p6<-ggplot(oplength,aes(x=value,fill=variable))+geom_density(alpha=0.4)+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p7<-ggplot(operonsize,aes(x=value,fill=variable))+geom_histogram(binwidth=3,alpha=0.5)+xlab("CDS per Operon")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))
ggsave("summary/images/assembly/utrlength.png",p1)
ggsave("summary/images/assembly/intrautrlength.png",p2)
ggsave("summary/images/assembly/igrlength.png",p3)
ggsave("summary/images/assembly/transcript_length.png",p4)
ggsave("summary/images/assembly/feature_length.png",p5)
ggsave("summary/images/assembly/operon_length.png",p6)
ggsave("summary/images/assembly/operon_size.png",p7)


#  PAIRED
putr<-rbind(data.frame(variable=rep("5`",length(putrfive$V2)),value=putrfive$V2),data.frame(variable=rep("3`",length(putrthree$V2)),value=putrthree$V2))
pintrautr<-rbind(data.frame(variable=rep("Ref. Intra-UTR",length(paredesintrautr$V2)),value=paredesintrautr$V2),data.frame(variable=rep("Intra-UTR",length(pintra$V2)),value=pintra$V2))

pigrs<-rbind(data.frame(variable=rep("Trinity",length(pigr$V3)),value=pigr$V3),data.frame(variable=rep("Ref. CDSes",length(prefigr$V3)),value=prefigr$V3),data.frame(variable=rep("Operons",length(paredesigr$V3)),value=paredesigr$V3))
plens<-rbind(data.frame(variable=rep("Novel",length(pnovel$V2)),value=pnovel$V2),data.frame(variable=rep("Standard",length(pstandard$V2)),value=pstandard$V2),data.frame(variable=rep("Total",length(ptotal)),value=ptotal))
ptotallens<-rbind(data.frame(variable=rep("Novel",length(pnovel$V2)),value=pnovel$V2),data.frame(variable=rep("Standard",length(pstandard$V2)),value=pstandard$V2),data.frame(variable=rep("Operon",length(poplen$V2)),value=poplen$V2),data.frame(variable=rep("Ref. CDS",length(prefcds$V2)),value=prefcds$V2),data.frame(variable=rep("Ref. Operon",length(paredesoplength$V1)),value=paredesoplength$V1),data.frame(variable=rep("Total",length(ptotal)),value=ptotal))

poplength<-rbind(data.frame(variable=rep("Ref. CDS",length(prefcds$V2)),value=prefcds$V2),data.frame(variable=rep("Ref. Operon",length(paredesoplength$V1[paredesopsize$V1 > 1])),value=paredesoplength$V1[paredesopsize$V1 > 1]),data.frame(variable=rep("Standard",length(pstandard$V2)),value=pstandard$V2),data.frame(variable=rep("Operon",length(poplen$V2)),value=poplen$V2))

poperonsize<-rbind(data.frame(variable=rep("Standard",length(popsize$V2)),value=popsize$V2),data.frame(variable=rep("Ref. Operons",length(paredesopsize$V1[paredesopsize$V1 > 1])),value=paredesopsize$V1[paredesopsize$V1 > 1]))


# Paired Graphing
# UTR length
p1<-ggplot(putr,aes(x=value,fill=variable))+geom_histogram(alpha=0.5,position="identity")+xlab("UTR length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Intra-UTR
p2.1<-ggplot(pintrautr,aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Intra-genic UTR length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+scale_x_log10(breaks=10**(0:3),labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='b')
p2.2<-ggplot(pintrautr[pintrautr$value > 0,],aes(x=value,fill=variable))+geom_histogram(alpha=0.9,position="dodge",binwidth=100)+scale_x_continuous(breaks=0:8*500,limits=c(0,4000))+xlab("Intra-genic UTR length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))
p2.3<-ggplot(pintrautr[pintrautr$value < 500,],aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Intra-genic UTR length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+scale_x_continuous(limits=c(-100,500))
# IGR length
p3<-ggplot(pigrs,aes(x=value,fill=variable))+geom_histogram(alpha=0.9,position="dodge")+xlab("IGR length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Novel, standard, total
p4<-ggplot(plens,aes(x=value,fill=variable))+geom_histogram(alpha=0.6,position="dodge")+xlab("Transcript length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Above plus CDS, ref operons, operons
p5<-ggplot(ptotallens,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("Feature length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Operon length
p6.1<-ggplot(poplength[(poplength$variable == "Ref. Operon" | poplength$variable == "Operon"),],aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Operon length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))
p6.2<-ggplot(poplength,aes(x=value,fill=variable))+geom_histogram(alpha=0.9,position="dodge")+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p6.3<-ggplot(poplength,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p6.4<-ggplot(poplength[(poplength$variable == "Ref. Operon" | poplength$variable == "Operon"),],aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Operon size
p7<-ggplot(poperonsize,aes(x=value,fill=variable))+geom_histogram(binwidth=2,alpha=0.5,position="identity")+xlab("CDS per Operon")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))
ggsave("summary/images/assembly/putrlength.png",p1)
ggsave("summary/images/assembly/pintrautrlength.1.png",p2.1)
ggsave("summary/images/assembly/pintrautrlength.2.png",p2.2)
ggsave("summary/images/assembly/pintrautrlength.3.png",p2.3)
ggsave("summary/images/assembly/pigrlength.png",p3)
ggsave("summary/images/assembly/ptranscript_length.png",p4)
ggsave("summary/images/assembly/pfeature_length.png",p5)
ggsave("summary/images/assembly/poperon_length.1.png",p6.1)
ggsave("summary/images/assembly/poperon_length.2.png",p6.2)
ggsave("summary/images/assembly/poperon_length.3.png",p6.3)
ggsave("summary/images/assembly/poperon_length.4.png",p6.4)
ggsave("summary/images/assembly/poperon_size.png",p7)







#  F I N A L   A N N O T A T I O N

futr<-rbind(data.frame(variable=rep("5`",length(futrfive$V2)),value=futrfive$V2),data.frame(variable=rep("3`",length(futrthree$V2)),value=futrthree$V2))
fintrautr<-rbind(data.frame(variable=rep("Ref. Intra-UTR",length(paredesintrautr$V2)),value=paredesintrautr$V2),data.frame(variable=rep("Intra-UTR",length(fintra$V2)),value=fintra$V2))

figrs<-rbind(data.frame(variable=rep("Trinity",length(figr$V3)),value=figr$V3),data.frame(variable=rep("Ref. CDSes",length(frefigr$V3)),value=frefigr$V3),data.frame(variable=rep("Operons",length(paredesigr$V3)),value=paredesigr$V3))
flens<-rbind(data.frame(variable=rep("Novel",length(fnovel$V2)),value=fnovel$V2),data.frame(variable=rep("Standard",length(fstandard$V2)),value=fstandard$V2),data.frame(variable=rep("Total",length(ftotal)),value=ftotal))
ftotallens<-rbind(data.frame(variable=rep("Novel",length(fnovel$V2)),value=fnovel$V2),data.frame(variable=rep("Standard",length(fstandard$V2)),value=fstandard$V2),data.frame(variable=rep("Operon",length(foplen$V2)),value=foplen$V2),data.frame(variable=rep("Ref. CDS",length(frefcds$V2)),value=frefcds$V2),data.frame(variable=rep("Ref. Operon",length(paredesoplength$V1)),value=paredesoplength$V1),data.frame(variable=rep("Total",length(ftotal)),value=ftotal))

foplength<-rbind(data.frame(variable=rep("Ref. CDS",length(frefcds$V2)),value=frefcds$V2),data.frame(variable=rep("Ref. Operon",length(paredesoplength$V1[paredesopsize$V1 > 1])),value=paredesoplength$V1[paredesopsize$V1 > 1]),data.frame(variable=rep("Standard",length(fstandard$V2)),value=fstandard$V2),data.frame(variable=rep("Operon",length(foplen$V2)),value=foplen$V2))

foperonsize<-rbind(data.frame(variable=rep("Standard",length(fopsize$V2)),value=fopsize$V2),data.frame(variable=rep("Ref. Operons",length(paredesopsize$V1[paredesopsize$V1 > 1])),value=paredesopsize$V1[paredesopsize$V1 > 1]))


# Paired Graphing
# UTR length
p1<-ggplot(futr,aes(x=value,fill=variable))+geom_histogram(alpha=0.5,position="identity")+xlab("UTR length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Intra-UTR
p2.1<-ggplot(fintrautr,aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Intra-genic UTR length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+scale_x_log10(breaks=10**(0:3),labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='b')
p2.2<-ggplot(fintrautr[fintrautr$value > 0,],aes(x=value,fill=variable))+geom_histogram(alpha=0.9,position="dodge",binwidth=100)+scale_x_continuous(breaks=0:8*500,limits=c(0,4000))+xlab("Intra-genic UTR length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+ggtitle("Intra-genic UTR > 0bp")
p2.3<-ggplot(fintrautr[fintrautr$value < 500,],aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Intra-genic UTR length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+scale_x_continuous(limits=c(-100,500))+ggtitle("Intra-genic UTR < 500bp")
# IGR length
p3.1<-ggplot(figrs,aes(x=value,fill=variable))+geom_histogram(alpha=0.9,position="dodge")+xlab("IGR length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p3.2<-ggplot(figrs,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("IGR length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Novel, standard, total
p4<-ggplot(flens,aes(x=value,fill=variable))+geom_histogram(alpha=0.6,position="dodge")+xlab("Transcript length")+scale_x_log10(limits=c(1,10**5),breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Above plus CDS, ref operons, operons
p5.1<-ggplot(ftotallens,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("Feature length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p5.2<-ggplot(ftotallens[ftotallens$value < 10000,],aes(x=value,fill=variable))+geom_density(alpha=0.4)+xlab("Feature length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+ggtitle("Transcript length < 10,000bp")
# Operon length
p6.1<-ggplot(foplength[(foplength$variable == "Ref. Operon" | foplength$variable == "Operon"),],aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Operon length")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))
p6.2<-ggplot(foplength,aes(x=value,fill=variable))+geom_histogram(alpha=0.9,position="dodge")+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p6.3<-ggplot(foplength,aes(x=value,fill=variable))+geom_density(alpha=0.5)+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
p6.4<-ggplot(foplength[(foplength$variable == "Ref. Operon" | foplength$variable == "Operon"),],aes(x=value,fill=variable))+geom_histogram(alpha=0.7,position="identity")+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
# Operon size
p7<-ggplot(foperonsize,aes(x=value,fill=variable))+geom_histogram(binwidth=2,alpha=0.5,position="identity")+xlab("CDS per Operon")+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))
ggsave("summary/images/assembly/futrlength.png",p1)
ggsave("summary/images/assembly/fintrautrlength.1.png",p2.1)
ggsave("summary/images/assembly/fintrautrlength.2.png",p2.2)
ggsave("summary/images/assembly/fintrautrlength.3.png",p2.3)
ggsave("summary/images/assembly/figrlength.1.png",p3.1)
ggsave("summary/images/assembly/figrlength.2.png",p3.2)
ggsave("summary/images/assembly/ftranscript_length.png",p4)
ggsave("summary/images/assembly/ffeature_length.1.png",p5.1)
ggsave("summary/images/assembly/ffeature_length.2.png",p5.2)
ggsave("summary/images/assembly/foperon_length.1.png",p6.1)
ggsave("summary/images/assembly/foperon_length.2.png",p6.2)
ggsave("summary/images/assembly/foperon_length.3.png",p6.3)
ggsave("summary/images/assembly/foperon_length.4.png",p6.4)
ggsave("summary/images/assembly/foperon_size.png",p7)








#  Final
p1<-ggplot(rbind(data.frame(variable=rep("Total",length(total)),value=total),data.frame(variable=rep("Paired-only",length(ptotal)),value=ptotal)),aes(x=value,fill=variable))+geom_density(alpha=0.4)+xlab("Operon length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))+theme(axis.text.y=element_text(colour="black"),axis.text.x=element_text(colour="black"))+annotation_logticks(base=10,sides='b')
ggsave("summary/images/pairedVtotal.png",p1)
