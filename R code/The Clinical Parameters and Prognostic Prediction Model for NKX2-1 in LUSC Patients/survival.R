#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("survival")
#install.packages("survminer")



library(limma)
library(survival)
library(survminer)

expFile="geneExp.txt"     
cliFile="time.txt"        
setwd("D:\\R code")     


rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
gene=colnames(rt)[1]


tumorData=rt[rt$Type=="Tumor",1,drop=F]
tumorData=as.matrix(tumorData)
rownames(tumorData)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(tumorData))
data=avereps(tumorData)


Type=ifelse(data[,gene]>median(data[,gene]), "High", "Low")
data=cbind(as.data.frame(data), Type)


cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli$futime=cli$futime/365


sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data)


outTab=cbind(ID=row.names(rt), rt)
outTab=outTab[,-ncol(outTab)]
write.table(outTab, file="expTime.txt", sep="\t", quote=F, row.names=F)


diff=survdiff(Surv(futime, fustat) ~ Type, data=rt)
pValue=1-pchisq(diff$chisq, df=1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=", sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ Type, data = rt)

		

surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           surv.median.line = "hv",
		           legend.title=gene,
		           legend.labs=c("High level", "Low level"),
		           xlab="Time(years)",
		           ylab="Overall survival",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=F,
		       	   risk.table.title="",
		           risk.table.col = "strata",
		           risk.table.height=.25)


pdf(file=paste0(gene, ".surv.pdf"), width=5.5, height=5, onefile=FALSE)
print(surPlot)
dev.off()




