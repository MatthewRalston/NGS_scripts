# Regularization.R
# Copyright 2014 Matt Ralston
# Updated: 7/27/2014
#
# This R script contains code to generate diagnostic and explanatory plots
# For the regularization process in DESeq2


# Variance vs mean trend (COV vs mean; Why do we regularize?)
mseq<-as.data.frame(mcols(deseq))
p1<-ggplot(mseq,aes(x=baseMean,y=sqrt(baseVar)/baseMean))+geom_point()+scale_x_log10(labels=trans_format('log10',math_format(10^.x)))+annotation_logticks(base=10,sides='b')+ylab(bquote(.("COV: ") ~ frac(sigma,mu)))+xlab(expression(mu))
ggsave("summary/images/COVvMean.png",p1)

# Dispersion plot (How do we regularize?)
p1<-ggplot(melt(as.data.frame(mcols(deseq))[,c(1,4,9,5)],id="baseMean"))+geom_point(aes(x=log10(baseMean),y=log10(value),colour=variable),alpha=0.7,size=1.5)+scale_colour_manual("Legend",values=c("black","#56B4E9","red"),labels=c("Dispersion",expression(paste("Max.",italic(' a posteriori'))),"Model fit"))+xlab(expression(paste(Log[10],"Regularized Expression")))+ylab(expression(paste(Log[10],"Dispersion")))+scale_x_continuous(labels=math_format(10^.x))+scale_y_continuous(breaks=-7:2,labels=math_format(10^.x))+annotation_logticks(sides='b')+theme(legend.justification=c(0,0),legend.position=c(0,0))+guides(colour=guide_legend(override.aes=list(size=5)))
ggsave("summary/images/Dispersion.png",p1)

# Residual plot (What is regularization?)
rawmodel<-lm(NS_A_15~NS_B_15,data=log10(rawcounts+1))
normmodel<-lm(NS_A_15~NS_B_15,data=log10(normcounts+1))
regmodel<-lm(NS_A_15~NS_B_15,data=(regcounts/log2(10)))
len<-length(residuals(rawmodel))
rawmodel<-data.frame(V1=rep("Raw counts",len),V2=fitted(rawmodel),V3=residuals(rawmodel))
normmodel<-data.frame(V1=rep("Norm. counts",len),V2=fitted(normmodel),V3=residuals(normmodel))
regmodel<-data.frame(V1=rep("Reg. counts",len),V2=fitted(regmodel),V3=residuals(regmodel))
resid<-rbind(rawmodel,normmodel,regmodel)
colnames(resid)<-c("Variable","Fitted","Residuals")
p1<-ggplot(resid)+geom_point(aes(x=Fitted,y=Residuals,colour=Variable),alpha=0.7,size=2)+scale_colour_manual("Legend",values=c("black","#56B4E9","red"),labels=c("Raw counts","Normalized counts", "Model fit"))+xlab("Gene Expression")+scale_y_continuous(breaks=-1:4,labels=math_format(10^.x))+ylab("Residuals")+scale_x_continuous(breaks=(-1:4),labels=math_format(10^.x))+annotation_logticks(base=10,sides='bl')
ggsave("summary/images/Residuals.png",p1)

# SD as a function of mean
library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(deseq))>0)
meanSdPlot(log2(counts(deseq,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
meanSdPlot(assay(rdeseq[notAllZero,]), ylim = c(0,2.5))
