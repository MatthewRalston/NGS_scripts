# Oldcompare.r
# 


                                        #        C I R C O S     S T U F F


x<-list(c(0.01,1.0),c(0.01,1.5),c(0.01,2.0),c(0.05,1.0),c(0.05,1.5),c(0.05,2.0))
# BuOH over time
buohgenes<-data.frame()
for (i in 1:length(buohtime)) {
    fc<-c()
    for (j in 1:length(x)) {
        a<-x[[j]][1]
        b<-x[[j]][2]
        tmp<-regcounts[rownames(subset(buohtime[[i]],subset=buohtime[[i]]$padj < a & abs(buohtime[[i]]$log2FoldChange) > b)),]
        fc<-c(fc,length(rownames(tmp)))
    }
    buohgenes<-rbind(buohgenes,fc)
}
colnames(buohgenes)<-c("1.0_0.01","1.5_0.01","2.0_0.01","1.0_0.05","1.5_0.05","2.0_0.05")
rownames(buohgenes)<-c("15","75","150","270")
write.table(buohgenes,file="summary/buoh-pairwise.txt",sep="\t",row.names=T,col.names=T,quote=F)

# Butyrate over time
bagenes<-data.frame()
for (i in 1:length(batime)) {
    fc<-c()
    for (j in 1:length(x)) {
        a<-x[[j]][1]
        b<-x[[j]][2]
        tmp<-regcounts[rownames(subset(batime[[i]],subset=batime[[i]]$padj < a & abs(batime[[i]]$log2FoldChange) > b)),]
        fc<-c(fc,length(rownames(tmp)))
    }
    bagenes<-rbind(bagenes,fc)
}

colnames(buohgenes)<-c("1.0_0.01","1.5_0.01","2.0_0.01","1.0_0.05","1.5_0.05","2.0_0.05")
rownames(buohgenes)<-c("15","75","150","270")
write.table(buohgenes,file="summary/butyrate-pairwise.txt",sep="\t",row.names=T,col.names=T,quote=F)

# Stress


buohfc<-c()
for (i in 1:length(x)) {
    a<-x[[i]][1]
    b<-x[[i]][2]
    buohfc<-c(buohfc,length(rownames(subset(stress[[1]],subset=stress[[1]]$padj < a & abs(stress[[1]]$log2FoldChange) > b))))
}
bafc<-c()
for (i in 1:length(x)) {
    a<-x[[i]][1]
    b<-x[[i]][2]
    bafc<-c(bafc,length(rownames(subset(stress[[2]],subset=stress[[2]]$padj < a & abs(stress[[2]]$log2FoldChange) > b))))
}
stressgenes<-rbind(buohfc,bafc)
colnames(stressgenes)<-c("1.0_0.01","1.5_0.01","2.0_0.01","1.0_0.05","1.5_0.05","2.0_0.05")
rownames(stressgenes)<-c("buoh","butyrate")
write.table(stressgenes,file="summary/general-pairwise.txt",sep="\t",row.names=T,col.names=T,quote=F)
