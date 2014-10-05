# ASSEMBLYSTATS.R
# Copyright 2014 Matt Ralston
# Updated: 10/5/14
#
# This is the R script to produce distributions of transcript length,etc for the assemblies.

utrfive<-read.table("Trinity_ref/stats/five.len",header=F)
utrthree<-read.table("Trinity_ref/stats/three.len",header=F)
igr<-read.table("Trinity_ref/stats/igr.len",header=F)
refigr<-read.table("Trinity_ref/stats/refigr.len",header=F)
opsize<-read.table("Trinity_ref/stats/operonsize.txt",header=F)
oplen<-read.table("Trinity_ref/stats/operon.len",header=F)
novel<-read.table("Trinity_ref/stats/novel.len",header=F)
standard<-read.table("Trinity_ref/stats/standard.len",header=F)
refcds<-read.table("Trinity_ref/stats/refcds.len",header=F)
total<-c(novel$V1,standard$V1)

utr<-rbind(data.frame(variable=rep("5`",length(utrfive$V1)),value=utrfive$V1),data.frame(variable=rep("3`",length(utrthree$V1)),value=utrthree$V1))

igrs<-rbind(data.frame(variable=rep("Trinity",length(igr$V1)),value=igr$V1),data.frame(variable=rep("Reference",length(refigr$V1)),value=refigr$V1))

lens<-rbind(data.frame(variable=rep("Novel",length(novel$V1)),value=novel$V1),data.frame(variable=rep("Standard",length(standard$V1)),value=standard$V1),data.frame(variable=rep("Operon",length(oplen$V1)),value=oplen$V1),data.frame(variable=rep("Reference",length(refcds$V1)),value=refcds$V1),data.frame(variable=rep("Total",length(total)),value=total))



p1<-ggplot(utr,aes(x=value,fill=variable))+geom_density(alpha=0.3)+xlab("UTR length")+scale_x_log10(breaks=10**(0:5),labels=trans_format("log10",math_format(10^.x)))
p2<-ggplot(igrs,aes(x=value,fill=variable))+geom_density(alpha=0.3)+xlab("IGR length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))
p3<-ggplot(lens,aes(x=value,fill=variable))+geom_density(alpha=0.7)+xlab("Feature length")+scale_x_log10(breaks=10**(0:6),labels=trans_format("log10",math_format(10^.x)))
p4<-ggplot(opsize,aes(x=V1))+geom_histogram(binwidth=3,colour="black",fill="white")+xlab("Operon Size")

ggsave("summary/images/utrlength.png",p1)
ggsave("summary/images/igrlength.png",p2)
ggsave("summary/images/transcript_length.png",p3)
ggsave("summary/images/operon_size.png",p4)
