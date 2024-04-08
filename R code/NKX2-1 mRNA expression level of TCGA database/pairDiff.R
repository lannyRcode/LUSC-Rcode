#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")



library(limma)
library(ggpubr)
expFile="geneExp.txt"      
setwd("D:\\R code")      


rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
geneName=colnames(rt)[1]


normalData=rt[rt$Type=="Normal",1,drop=F]
normalData=as.matrix(normalData)
rownames(normalData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(normalData))
normalData=avereps(normalData)
tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
tumorData=avereps(tumorData)
sameSample=intersect(row.names(normalData), row.names(tumorData))
data=cbind(normalData[sameSample,,drop=F], tumorData[sameSample,,drop=F])
colnames(data)=c("Normal", "Tumor")
data=as.data.frame(data)


pdf(file="pairDiff.pdf", width=5, height=4.5)
ggpaired(data, cond1="Normal", cond2="Tumor", fill="condition",
	xlab="", ylab=paste0(geneName, " expression"),
	legend.title="Type",
	palette=c("blue","red"))+
    stat_compare_means(paired = TRUE, symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif",label.x = 1.35)
dev.off()



