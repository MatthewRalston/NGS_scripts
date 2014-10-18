

assembly<-read.table("Trinity_ref/Trinity_filtered.gtf",header=F)
alength<-assembly[,5]-assembly[,4]
reference<-read.table("reference/CAC_proteins.gtf",header=F)
rlength<-reference[,5]-reference[,4]

rlength<-data.frame(variable=rep("ref",length(rlength)),value=rlength)
alength<-data.frame(variable=rep("trin",length(alength)),value=alength)
lens<-rbind(alength,rlength)

p1<-ggplot(lens,aes(x=value,fill=variable))+geom_density(alpha=.3)+scale_x_log10(labels=trans_format("log10",math_format(10^.x)))+annotation_logticks(base=10,sides="b")+xlab("Transcript Size")
ggsave("summary/images/assembly_summary.png",p1)
