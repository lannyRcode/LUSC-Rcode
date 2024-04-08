#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


library(limma)
library(ggpubr)

expFile="geneExp.txt"       
cliFile="clinical.txt"      
setwd("D:\\R code")     


rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]


tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)

cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))


samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)


for(clinical in colnames(rt[,2:ncol(rt)])){
	data=rt[c(gene, clinical)]
	colnames(data)=c(gene, "clinical")
	data=data[(data[,"clinical"]!="unknow"),]
	group=levels(factor(data$clinical))
	data$clinical=factor(data$clinical, levels=group)
	comp=combn(group,2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
		          xlab=clinical,
		          ylab=paste(gene, " expression"),
		          legend.title=clinical)+ 
	    stat_compare_means(comparisons = my_comparisons)
	pdf(file=paste0("clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
	print(boxplot)
	dev.off()
}




