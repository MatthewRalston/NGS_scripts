# AlignmentStats.R
# Copyright 2014 Matt Ralston
# Updated: 7/26/14
#
# This R script contains code to generate summary plots for the alignment
# of Fastq reads to the C_ac genome.

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


