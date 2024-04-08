#source("https://bioconductor.org/biocLite.R")
#biocLite("sva")
#biocLite("limma")

library(sva)
library(limma)
setwd("D:\\R code")
rt=read.table("merge.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
batchType=c(rep(1,129),rep(2,18))
modType=c(rep("normal",60),rep("tumor",69),rep("normal",9),rep("tumor",9))
mod = model.matrix(~as.factor(modType))
outTab=ComBat(data, batchType, mod, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)



