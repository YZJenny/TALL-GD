####################
## 3. 看gene的表达(cutoff由surv_cutpoint设定)在survial(CLIN/OS)的关系
####################
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(ggplot2)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

EXP <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/inputExp.txt',header = T,sep='\t')
EFS <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/EFS_subtype.txt',header = T,row.names = 1,sep='\t')
OS <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/OS_subtype.txt',header = T,row.names = 1,sep='\t')

Type <- c('EFS','OS')
pvalue_df <- c()

for(j in 1:length(Type)){
  type=Type[j]
  if(type=='EFS'){
    CLIN=EFS
  }else if(type=='OS'){
    CLIN=OS
  }
  CLIN$surtype <- NULL
  
  EXP_tmp <- as.data.frame(t(EXP))
  EXP_tmp <- EXP_tmp[,sapply(EXP_tmp, function(x)
    ifelse(sum(x != 0, na.rm = TRUE) > 30, TRUE, FALSE)) ] #删除列中包含>30个0的基因，否则会报错
  EXP_tmp <- tibble::rownames_to_column(EXP_tmp,'ID')
  CLIN_tmp <- tibble::rownames_to_column(CLIN,'ID')
  data <- merge(CLIN_tmp,EXP_tmp)
  
  Genes=colnames(data)[4:ncol(data)]
  res.cut <- surv_cutpoint(data,time = "time", event = "status",variables = Genes)
  cut.df <- as.data.frame(summary(res.cut))
  cut.df <- tibble::rownames_to_column(cut.df,'gene')
  
  res.cat <- surv_categorize(res.cut)
  rownames(res.cat) <- data$ID
  write.table(cut.df,
              paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/Cutpoint_',type,'.txt',sep=''),
              row.names = F,col.names = T,quote = F,sep='\t')
  
  write.table(res.cat,
              paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/Cat_',type,'.txt',sep=''),
              row.names = T,col.names = T,quote = F,sep='\t')
  
  clin_info=res.cat[,1:2]
  group_info <- res.cat[,3:ncol(res.cat)]
  pvalue <- c()
  for(k in 1:length(Genes)){
    print(k)
    surtype=group_info[,k]
    clin_info$surtype <- surtype
    fit <- survfit(Surv(time =time, event =status) ~ surtype,data=clin_info)
    pvalue[k] <- surv_pvalue(fit)$pval
  }
  names(pvalue) <- Genes
  print(length(pvalue[pvalue <0.05]))
  pvalue_df <- rbind(pvalue_df,round(pvalue,3))
}

rownames(pvalue_df) <- Type
pvalue_df <- as.data.frame(t(pvalue_df))
pvalue_df <- tibble::rownames_to_column(pvalue_df,'gene')

write.table(pvalue_df,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/Pvalue.txt',
            col.names = T,row.names = F,quote = F,sep='\t')

# pvalue_df  <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/Pvalue.txt',
#                          header  = T,sep='\t')
head(pvalue_df)
Surv.EFS <- pvalue_df$gene[pvalue_df$EFS<=0.05]
Surv.OS <- pvalue_df$gene[pvalue_df$OS<=0.05]
geneLst <- list(EFS=Surv.EFS,OS=Surv.OS)
saveRDS(geneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/sigSurv.RDS')

geneSet <- union(Surv.OS,Surv.EFS)
saveRDS(geneSet,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/GeneLst.RDS')

##根据uni-cox的HR 判断方向
# Unires.EFS <- as.data.frame(readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/Unires_EFS.RDS'))
# Unires.EFS <- Unires.EFS[-1,]
# Unires.EFS <- Unires.EFS$Gene[Unires.EFS$HR>1]
# 
# Unires.OS <- as.data.frame(readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/UniCox/Unires_OS.RDS'))
# Unires.OS <- Unires.OS[-1,]
# Unires.OS <- Unires.OS$Gene[Unires.OS$HR>1]
# 
# print(length(Surv.EFS));print(length(Surv.OS))
# print(length(Unires.EFS));print(length(Unires.OS))
# EFS <- intersect(Surv.EFS,Unires.EFS);length(EFS)
# OS <- intersect(Surv.OS,Unires.OS);length(OS)

