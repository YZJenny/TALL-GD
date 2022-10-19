######
###  UniCOX Analysis for gene
######
rm(list=ls())
library(ggplot2)
library(pheatmap)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
library(rmda)
library(DCA)
get_Unires <- function(Cox_df){
  covariates <- colnames(Cox_df)[4:ncol(Cox_df)]
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(time, status)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = Cox_df)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$sctest["pvalue"], digits=2)
                           sc.test<-signif(x$sctest["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR.confint <- paste0(HR, " (", 
                                                HR.confint.lower, "-", HR.confint.upper, ")")
                           # res<-c(beta, HR.confint, sc.test, p.value)
                           # names(res)<-c("beta", "HR (95% CI for HR)", "sc.test", "p.value")
                           res<-c(HR, HR.confint.lower, HR.confint.upper, HR.confint, p.value)
                           names(res)<-c( 'HR','lower','upper',"HR (95% CI for HR)", "p.value")
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res,stringsAsFactors = F)
  # res$Gene <- rownames(res)
  # res$adjustP <- p.adjust(res$p.value,method = 'fdr')
  # res <- as.matrix(res)
  # res <- rbind(colnames(res),res)
  return(res)
}


## survival file
Type='EFS'
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

Cox_df <- merge(Surv,sub.EXP,by='sampleID')
dim(Cox_df)

Unires <- get_Unires(Cox_df)
sig.Unires <- Unires[Unires$p.value<0.05,]
dim(sig.Unires)

system('mkdir -p /local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/')
saveRDS(Unires,paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/Unires_',Type,'.RDS',sep=''))
saveRDS(sig.Unires,paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/sigUnires_',Type,'.RDS',sep=''))


favorG=rownames(sig.Unires)[sig.Unires$HR<1]
unfavorG=rownames(sig.Unires)[sig.Unires$HR>1]
progG.lst=list(favorG=favorG,unfavorG=unfavorG)
saveRDS(progG.lst,paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/progGeneLst_',Type,'.RDS',sep=''))
