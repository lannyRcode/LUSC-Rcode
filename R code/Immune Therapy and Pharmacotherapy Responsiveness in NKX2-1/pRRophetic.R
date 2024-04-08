#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("limma", "car", "ridge", "preprocessCore", "genefilter", "sva"))

#install.packages("ggplot2")
#install.packages("ggpubr")


library(limma)
library(ggpubr)
library(pRRophetic)
library(ggplot2)
set.seed(12345)

pFilter=0.001            
gene="NKX2-1"             
expFile="symbol.txt"     
setwd("D:\\R code")     

data(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
allDrugs=unique(drugData2016$Drug.name)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0.5,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2","1",group)
data=data[,group==0]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=t(avereps(data))

geneExp=as.data.frame(t(data[gene,,drop=F]))
geneExp$Type=ifelse(geneExp[,gene]>median(geneExp[,gene]), "High", "Low")

for(drug in allDrugs){
	possibleError=tryCatch(
    	{senstivity=pRRopheticPredict(data, drug, selection=1, dataset = "cgp2016")},
    	error=function(e) e)
    if(inherits(possibleError, "error")){next}
	senstivity=senstivity[senstivity!="NaN"]
	senstivity[senstivity>quantile(senstivity,0.99)]=quantile(senstivity,0.99)
	
	sameSample=intersect(row.names(geneExp), names(senstivity))
	geneExp=geneExp[sameSample, "Type",drop=F]
	senstivity=senstivity[sameSample]
	rt=cbind(geneExp, senstivity)
	
	rt$Type=factor(rt$Type, levels=c("Low", "High"))
	type=levels(factor(rt[,"Type"]))
	comp=combn(type, 2)
	my_comparisons=list()
	for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
	
	test=wilcox.test(senstivity~Type, data=rt)
	diffPvalue=test$p.value
	if(diffPvalue<pFilter){
		boxplot=ggboxplot(rt, x="Type", y="senstivity", fill="Type",
					      xlab=gene,
					      ylab=paste0(drug, " senstivity (IC50)"),
					      legend.title=gene,
					      palette=c("#0066FF","#FF0000")
					     )+ 
			stat_compare_means(comparisons=my_comparisons)
		pdf(file=paste0("durgSenstivity.", drug, ".pdf"), width=5, height=4.5)
		print(boxplot)
		dev.off()
	}
}



