#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("ggExtra")


#ÒýÓÃ°ü
library(limma)
library(ggplot2)
library(ggpubr)
library(ggExtra)

gene="NKX2-1"          
corFilter=0.6             
pFilter=0.001             
expFile="symbol.txt"      
setwd("D:\\R code")      


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>1,]

group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0]
data=log2(data+1)


x=as.numeric(data[gene,])
outTab=data.frame()
for(j in rownames(data)){
	if(gene==j){next}
    y=as.numeric(data[j,])
	corT=cor.test(x, y, method = 'pearson')
	cor=corT$estimate
	pvalue=corT$p.value
	outTab=rbind(outTab, cbind(Query=gene, Gene=j, cor, pvalue))
	
	if((abs(cor)>corFilter) & (pvalue<pFilter)){
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
			xlab(paste0(gene, " expression"))+ ylab(paste0(j, " expression"))+
			geom_point()+ geom_smooth(method="lm", formula=y~x) + theme_bw()+
			stat_cor(method = 'pearson', aes(x =x, y =y))
		pdf(file=paste0("cor.", j, ".pdf"), width=5, height=4.6)
		print(p1)
		dev.off()
	}
}


write.table(file="corResult.txt", outTab, sep="\t", quote=F, row.names=F)
outTab=outTab[abs(as.numeric(outTab$cor))>corFilter & as.numeric(outTab$pvalue)<pFilter,]
write.table(file="corSig.txt", outTab, sep="\t", quote=F, row.names=F)




