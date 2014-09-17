# PCA.R
# Copyright 2014 Matt Ralston
# Updated: 4/27/2014
#
# This script contains code for principal component analysis of expression data
require('rgl')
require('gdata')

# PCA of whole data
pca<-princomp(regcounts,scale=TRUE)
cormat<-cor(regcounts)
covmat<-cov(regcounts)
ba<-colnames(regcounts)[c(2,8,5,10,6,1,7,3,9,4)]
buoh<-colnames(regcounts)[c(12,18,15,20,16,11,17,13,19,14)]
ns<-colnames(regcounts)[c(22,28,25,30,26,21,27,23,29,24)]
bacol<-c("#99FFCC","#99FFCC","#33FF66","#33FF66","#33FF66","#33FF33","#33FF33","#33FF00","#33FF00","#33FF00")
buohcol<-c("#6600CC","#6600CC","#CC00CC","#CC00CC","#CC00CC","#FF3399","#FF3399","#FF0066","#FF0066","#FF0066")
nscol<-c("#00FFFF","#00FFFF","#0099FF","#0099FF","#0099FF","#3366CC","#3366CC","#3333CC","#3333CC","#3333CC")
regcolors<-c(rep("green",10),c(rep("red",10)),rep("blue",10))
corscores<-data.frame(samples=colnames(regcounts),prcomp(regcounts,scale=TRUE,center=TRUE)$rotation)[c(2,8,5,10,6,1,7,3,9,4,12,18,15,20,16,11,17,13,19,14,22,28,25,30,26,21,27,23,29,24),]
covscores<-data.frame(samples=colnames(regcounts),prcomp(covmat,scale=TRUE)$rotation)[c(2,8,5,10,6,1,7,3,9,4,12,18,15,20,16,11,17,13,19,14,22,28,25,30,26,21,27,23,29,24),]

# Basic pca
plot(pca)
biplot(pca)
# Heatmaps of PC1 and PC2
p1<-ggplot(melt(cormat))+geom_tile(aes(Var1,Var2,fill=value))+ggtitle("Correlation")
p2<-ggplot(melt(covmat))+geom_tile(aes(Var1,Var2,fill=value))+ggtitle("Covariance")
p<-arrangeGrob(p1,p2,ncol=2)
ggsave("summary/images/pca/heatmaps.png")
# PCA of correlation matrix
corvar=prcomp(cormat,scale=TRUE)$sdev**2
cornames=paste("PC",(1:length(corvar)),sep='')
cornames=factor(cornames,levels=cornames)
corscree<-data.frame(variable=cornames,value=corvar)
p1<-ggplot(corscree)+geom_bar(aes(variable,value),stat='identity')+ylab("Variance")+theme(axis.title.x=element_blank())+ggtitle("Principal Components of Correlation Matrix")
ggsave("summary/images/pca/cor_scree.png",p1)
# PCA of covariance matrix
covvar=prcomp(covmat,scale=TRUE)$sdev**2
covnames=paste("PC",(1:length(covvar)),sep='')
covnames=factor(covnames,levels=covnames)
covscree<-data.frame(variable=covnames,value=covvar)
p2<-ggplot(covscree)+geom_bar(aes(variable,value),stat='identity')+ylab("Variance")+theme(axis.title.x=element_blank())+ggtitle("Princpal Components of Covariance Matrix")
ggsave("summary/images/pca/cov_scree.png",p2)


# SIMPLE PC plots
# PC1 v PC2
#    correlation
#corscores$colors<-c(bacol,buohcol,nscol)
corscores$colors<-regcolors
p1<-ggplot(corscores)+geom_point(aes(x=PC1,y=PC2,colour=factor(samples,levels=samples)),stat='identity',size=5)+scale_colour_manual("Legend",values=corscores$colors)
ggsave("summary/images/pca/cor_PC1vPC2.png",p1)
#     covariance
#covscores$colors<-c(bacol,buohcol,nscol)
covscores$colors<-regcolors
p2<-ggplot(covscores)+geom_point(aes(x=PC1,y=PC2,colour=factor(samples,levels=samples)),stat='identity',size=5)+scale_colour_manual("Legend",values=covscores$colors)
ggsave("summary/images/pca/cov_PCvPC2.png",p2)
# PC2 v PC3
corscores$colors<-regcolors
p1<-ggplot(corscores)+geom_point(aes(x=PC2,y=PC3,colour=factor(samples,levels=samples)),stat='identity',size=5)+scale_colour_manual("Legend",values=corscores$colors)
ggsave("summary/images/pca/cor_PC2vPC3.png",p1)
#     covariance
#covscores$colors<-c(bacol,buohcol,nscol)
covscores$colors<-regcolors
p2<-ggplot(covscores)+geom_point(aes(x=PC1,y=PC3,colour=factor(samples,levels=samples)),stat='identity',size=5)+scale_colour_manual("Legend",values=covscores$colors)
ggsave("summary/images/pca/cov_PC2vPC3.png",p2)
# PC1 v PC3
corscores$colors<-regcolors
p1<-ggplot(corscores)+geom_point(aes(x=PC1,y=PC3,colour=factor(samples,levels=samples)),stat='identity',size=5)+scale_colour_manual("Legend",values=corscores$colors)
ggsave("summary/images/pca/cor_PC1vPC3.png",p1)
#     covariance
#covscores$colors<-c(bacol,buohcol,nscol)
covscores$colors<-regcolors
p2<-ggplot(covscores)+geom_point(aes(x=PC1,y=PC3,colour=factor(samples,levels=samples)),stat='identity',size=5)+scale_colour_manual("Legend",values=covscores$colors)
ggsave("summary/images/pca/cov_PC1vPC3.png",p2)





# 3D, interactive PCA
#    PCs 1, 2, and 3
scores<-covscores

x<-scores[,2:4]
zeroes<-rep(0,length(x$PC1))
x<-interleave(data.frame(PC1=zeroes,PC2=zeroes,PC3=zeroes),x)
rgl.open()
rgl.bg(color="white")
z<-x[1:20,]
rgl.lines(z[,1],z[,2],z[,3],color='green')
z<-x[21:40,]
rgl.lines(z[,1],z[,2],z[,3],color='red')
z<-x[41:60,]
rgl.lines(z[,1],z[,2],z[,3],color='blue')
axes3d(c('x','y','z'),col='black',labels=FALSE)
title3d('','','PC1','PC2','PC3')
n<-c(paste(rep("00",9),1:9,sep=''),paste(rep("0",89),10:99,sep=''),as.character(100:360))
for (i in 1:length(n)) {
    rgl.viewpoint(i,0,interactive=F)
    rgl.snapshot(paste("summary/pca/panoramic/file",n[i],".png",sep=''))
}
system('convert -delay 1 summary/pca/panoramic/*.png summary/images/pca/pca123.gif')
writeWebGL(dir="summary/pca/webgl/PC123")
#    PCs 1,2,4
x<-scores[,c(2,3,5)]
zeroes<-rep(0,length(x$PC1))
x<-interleave(data.frame(PC1=zeroes,PC2=zeroes,PC4=zeroes),x)
rgl.open()
rgl.bg(color="white")
z<-x[1:20,]
rgl.lines(z[,1],z[,2],z[,3],color='green')
z<-x[21:40,]
rgl.lines(z[,1],z[,2],z[,3],color='red')
z<-x[41:60,]
rgl.lines(z[,1],z[,2],z[,3],color='blue')
axes3d(c('x','y','z'),col='black',labels=FALSE)
title3d('','','PC1','PC2','PC4')
for (i in 1:length(n)) {
    rgl.viewpoint(i,0,interactive=F)
    rgl.snapshot(paste("summary/pca/panoramic/file",n[i],".png",sep=''))
}
system('convert -delay 1 summary/pca/panoramic/*.png summary/images/pca/pca124.gif')
writeWebGL(dir="summary/pca/webgl/PC124")

#    PCs 1,3,4
x<-scores[,c(2,4,5)]
zeroes<-rep(0,length(x$PC1))
x<-interleave(data.frame(PC1=zeroes,PC3=zeroes,PC4=zeroes),x)
rgl.open()
rgl.bg(color="white")
z<-x[1:20,]
rgl.lines(z[,1],z[,2],z[,3],color='green')
z<-x[21:40,]
rgl.lines(z[,1],z[,2],z[,3],color='red')
z<-x[41:60,]
rgl.lines(z[,1],z[,2],z[,3],color='blue')
axes3d(c('x','y','z'),col='black',labels=FALSE)
title3d('','','PC1','PC3','PC4')
for (i in 1:length(n)) {
    rgl.viewpoint(i,0,interactive=F)
    rgl.snapshot(paste("summary/pca/panoramic/file",n[i],".png",sep=''))
}
system('convert -delay 1 summary/pca/panoramic/*.png summary/images/pca/pca134.gif')
writeWebGL(dir="summary/pca/webgl/PC134")


#    PCs 2,3,4
x<-scores[,3:5]
zeroes<-rep(0,length(x$PC2))
x<-interleave(data.frame(PC2=zeroes,PC3=zeroes,PC4=zeroes),x)
rgl.open()
rgl.bg(color="white")
z<-x[1:20,]
rgl.lines(z[,1],z[,2],z[,3],color='green')
z<-x[21:40,]
rgl.lines(z[,1],z[,2],z[,3],color='red')
z<-x[41:60,]
rgl.lines(z[,1],z[,2],z[,3],color='blue')
axes3d(c('x','y','z'),col='black',labels=FALSE)
title3d('','','PC2','PC3','PC4')

n<-c(paste(rep("00",9),1:9,sep=''),paste(rep("0",89),10:99,sep=''),as.character(100:360))
for (i in 1:length(n)) {
    rgl.viewpoint(i,0,interactive=F)
    rgl.snapshot(paste("summary/pca/panoramic/file",n[i],".png",sep=''))
}
system('convert -delay 1 summary/pca/panoramic/*.png summary/images/pca/pca234.gif')
writeWebGL(dir="summary/pca/webgl/PC234")




