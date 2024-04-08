#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

expFile="geneExp.txt"     
tmbFile="TMB.txt"         
setwd("D:\\R code")     

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]

tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3\\-\\4", rownames(tumorData))
data=avereps(tumorData)

tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)

sameSample=intersect(row.names(data), row.names(tmb))
data=data[sameSample,,drop=F]
tmb=tmb[sameSample,,drop=F]
rt=cbind(data, tmb)

x=as.numeric(rt[,gene])
y=log2(as.numeric(rt[,"TMB"])+1)
df1=as.data.frame(cbind(x,y))
corT=cor.test(x, y, method="spearman")
p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0(gene, " expression"))+ylab("Tumor mutation burden")+
			geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
			stat_cor(method = 'spearman', aes(x =x, y =y))
p2=ggMarginal(p1, type = "density", xparams = list(fill = "orange"),yparams = list(fill = "blue"))

pdf(file="cor.pdf",width=5,height=5)
print(p2)
dev.off()



