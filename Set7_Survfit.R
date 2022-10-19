######
###  Survfit Analysis for gene
######
rm(list=ls())
library(ggplot2)
library(pheatmap)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
library(rmda)


## survival file
Type='OS'
Surv <- read.table(paste('/local/yanzijun/CRU/TALL/data/NG_2017/',Type,'_TAL.txt',sep=''),header = T,row.names = 1,sep='\t')
Surv <- tibble::rownames_to_column(Surv,'sampleID')
Surv <- Surv[,1:3]
head(Surv)


## gene exp file
EXP <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/inputExp.txt',header = T,sep='\t')
if(max(EXP)>50){
  EXP <- log(EXP+1,2)
}
max(EXP)
min(EXP)
EXP <- dplyr::select(EXP,Surv$sampleID)
dim(EXP)
print(all(colnames(EXP)==Surv$sampleID))
EXP <- as.data.frame(t(EXP))
EXP$sampleID <- rownames(EXP)

sub.EXP <- EXP[,-c(grep('-',colnames(EXP)),
                   grep('@',colnames(EXP)),
                   grep('\\.',colnames(EXP)),
                   grep('7A5',colnames(EXP)),
                   grep('NOP5',colnames(EXP)))]
dim(sub.EXP)

print(all(sub.EXP$sampleID==Surv$sampleID))

## 表达值变成二分类
tmp <- function(x){
  survtype <- rep('Low',length(x))
  survtype[x>quantile(x,0.5)] <- 'High'
  return(survtype)
}
tmpp <- apply(sub.EXP[,-ncol(sub.EXP)],2,tmp)


Cox_df <- cbind(Surv,tmpp)
dim(Cox_df)

## 计算p value
pvalue.lst <- c()
for(gene in colnames(Cox_df)[4:ncol(Cox_df)]){
  print(gene)
  surv.df <- Cox_df[,c('status','time',gene)]
  colnames(surv.df)[3] <- 'surtype'
  fit <- survfit(Surv(time =time, event =status) ~ surtype,data=surv.df)
  pvalue <- surv_pvalue(fit)$pval
  pvalue.lst <- c(pvalue.lst,pvalue)
}
names(pvalue.lst) <- colnames(Cox_df)[4:ncol(Cox_df)]


sig.pvalue <- pvalue.lst[pvalue.lst<0.05]
length(sig.pvalue)

system('mkdir -p /local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/')
saveRDS(pvalue.lst,paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/Survfit_',Type,'.RDS',sep=''))
saveRDS(sig.pvalue,paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/sigSurvfit_',Type,'.RDS',sep=''))
