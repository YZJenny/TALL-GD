############################################################
### Set1:DEG, get数据集的DEG和met_express的输入矩阵
############################################################
rm(list=ls())
library(readxl)
library(limma)

##############################
### GSE26713
##############################
meta <- read.csv('/local/yanzijun/CRU/TALL/data/GSE26713/metaInfo.txt',sep='\t')
data <- read.table('/local/yanzijun/CRU/TALL/data/GSE26713/exprSet.txt',sep='\t')
data <- dplyr::select(data,meta$geo_accession)
data <- data[order(rownames(data),decreasing = F),] ##基因按顺序排好，为何会有纯数字的基因，删除

print(all(colnames(data)==meta$geo_accession))
label <- meta$cytogenetics.ch1
label[label!='BM']='tumor'
label[label=='BM']='normal'
Group <- data.frame(ID=colnames(data),ID_REF=label)

write.table(data,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/Gene_Expression_Value.txt',quote = F,sep='\t',row.names = T,col.names = T)
write.table(Group,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/DEG_design.txt',quote = F,sep='\t',row.names = F,col.names = F)


## 执行limma的DEG
label <- meta$cytogenetics.ch1
label[label!='BM']='TALL'
label[label=='BM']='Normal'
EXP=log(data+1,2)
group_list <- factor(label,levels=c('Normal','TALL'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(TALL-Normal, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene')  
print(output[output$gene %in% c('TAL1','S1PR3','ARRB1'),])
write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/Limma_output.txt',sep='\t',row.names = F,quote = F)

## DEGlst
up.Genes <- output$gene[output$logFC > 0 & output$adj.P.Val < 0.05]
dn.Genes <- output$gene[output$logFC < 0 & output$adj.P.Val < 0.05]
print(length(up.Genes))
print(length(dn.Genes))
GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/GeneLst.RDS')

## DEG for met_express
#output <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/Limma_output.txt',sep='\t',header = T)
sig.output <- output[abs(output$logFC)> log(1,2) & output$adj.P.Val<0.05,]
write.table(sig.output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/DEG_P0.05_Lfc0.txt',sep='\t',row.names = F,quote = F)



##############################
### GSE48558
##############################
meta <- read.csv('/local/yanzijun/CRU/TALL/data/GSE48558/metaInfo.txt',sep='\t')
data <- read.table('/local/yanzijun/CRU/TALL/data/GSE48558/exprSet.txt',sep='\t')

## 挑选T-ALL病人
meta <- meta[meta$source_name_ch1 %in% c('T ALL Patient','T Cells'),]
data <- dplyr::select(data,rownames(meta))
dim(data)

data <- data[order(rownames(data),decreasing = F),] ##基因按顺序排好，为何会有纯数字的基因，删除

print(all(colnames(data)==rownames(meta)))

label <- meta$source_name_ch1
label[label=='T Cells']='normal'
label[label=='T ALL Patient']='tumor'
Group <- data.frame(ID=colnames(data),ID_REF=label)

write.table(data,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/Gene_Expression_Value.txt',quote = F,sep='\t',row.names = T,col.names = T)
write.table(Group,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/DEG_design.txt',quote = F,sep='\t',row.names = F,col.names = F)


## 执行limma的DEG
print(max(data))

## 判断是否log转换
if(max(data) > 50){
  EXP=log(data+1,2)
}else{
  EXP=data
}
group_list <- factor(label,levels=c('normal','tumor'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene')  
print(output[output$gene %in% c('TAL1','S1PR3','ARRB1'),])
write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/Limma_output.txt',sep='\t',row.names = F,quote = F)

## DEGlst
up.Genes <- output$gene[output$logFC > 0 & output$adj.P.Val < 0.05]
dn.Genes <- output$gene[output$logFC < 0 & output$adj.P.Val < 0.05]
print(length(up.Genes))
print(length(dn.Genes))
GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/GeneLst.RDS')

## DEG for met_express
#output <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/Limma_output.txt',sep='\t',header = T)
sig.output <- output[abs(output$logFC)> log(1,2) & output$adj.P.Val<0.05,]
write.table(sig.output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/DEG_P0.05_Lfc0.txt',sep='\t',row.names = F,quote = F)



##############################
### GSE66638(8 TALL/2 thymus)和GSE50999(43 TALL)
##############################
rm(list=ls())
library(readxl)
library(limma)
#### 合并两套数据
data1=read.csv('CRU/TALL/data/GSE50999/exprSet.txt',header = T,sep='\t')
label1 <- rep('tumor',ncol(data1))
data1$gene=rownames(data1)
data2=read.csv('CRU/TALL/data/GSE66638/exprSet.txt',header = T,sep='\t')
meta2 <- read.csv('/local/yanzijun/CRU/TALL/data/GSE66638/metaInfo.txt',sep='\t')
dim(meta2)
data2 <- dplyr::select(data2,rownames(meta2))
data2$gene=rownames(data2)

data=merge(data1,data2,by='gene')
rownames(data)=data$gene;data$gene=NULL
dim(data)

data <- data[order(rownames(data),decreasing = F),] ##基因按顺序排好，为何会有纯数字的基因，删除

label2 <- meta2$title
label2[grep('^T-ALL',label2)]='tumor'
label2[grep('^Thymocytes',label2)]='normal'
print(table(label2))
label=c(label1,label2)
print(table(label))
Group <- data.frame(ID=colnames(data),ID_REF=label)

write.table(data,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE66638_GSE50999/Gene_Expression_Value.txt',quote = F,sep='\t',row.names = T,col.names = T)
write.table(Group,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE66638_GSE50999/DEG_design.txt',quote = F,sep='\t',row.names = F,col.names = F)


## 执行limma的DEG
print(max(data))

## 判断是否log转换
if(max(data) > 50){
  EXP=log(data+1,2)
}else{
  EXP=data
}
group_list <- factor(label,levels=c('normal','tumor'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene')  
print(output[output$gene %in% c('TAL1','S1PR3','ARRB1'),])
write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE66638_GSE50999/Limma_output.txt',sep='\t',row.names = F,quote = F)

## DEGlst
up.Genes <- output$gene[output$logFC > 0 & output$adj.P.Val < 0.05]
dn.Genes <- output$gene[output$logFC < 0 & output$adj.P.Val < 0.05]
print(length(up.Genes))
print(length(dn.Genes))
GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE66638_GSE50999/GeneLst.RDS')

## DEG for met_express
#output <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE66638/Limma_output.txt',sep='\t',header = T)
sig.output <- output[abs(output$logFC)> log(1,2) & output$adj.P.Val<0.05,]
write.table(sig.output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE66638_GSE50999/DEG_P0.05_Lfc0.txt',sep='\t',row.names = F,quote = F)


##############################
### GSE41621(17 TALL/2 Tcell)
##############################
rm(list=ls())
library(readxl)
library(limma)

meta <- read.csv('/local/yanzijun/CRU/TALL/data/GSE41621/metaInfo.txt',sep='\t')
data <- read.table('/local/yanzijun/CRU/TALL/data/GSE41621/exprSet.txt',sep='\t')
dim(meta)
dim(data)

data <- dplyr::select(data,rownames(meta))
dim(data)

data <- data[order(rownames(data),decreasing = F),] ##基因按顺序排好，为何会有纯数字的基因，删除

print(all(colnames(data)==rownames(meta)))

label <- meta$title
label[grep('^T-ALL',label)]='tumor'
label[grep('^Control',label)]='normal'
print(table(label))
Group <- data.frame(ID=colnames(data),ID_REF=label)

write.table(data,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE41621/Gene_Expression_Value.txt',quote = F,sep='\t',row.names = T,col.names = T)
write.table(Group,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE41621/DEG_design.txt',quote = F,sep='\t',row.names = F,col.names = F)


## 执行limma的DEG
print(max(data))

## 判断是否log转换
if(max(data) > 50){
  EXP=log(data+1,2)
}else{
  EXP=data
}
group_list <- factor(label,levels=c('normal','tumor'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene')  
print(output[output$gene %in% c('TAL1','S1PR3','ARRB1'),])
write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE41621/Limma_output.txt',sep='\t',row.names = F,quote = F)

## DEGlst
up.Genes <- output$gene[output$logFC > 0 & output$adj.P.Val < 0.05]
dn.Genes <- output$gene[output$logFC < 0 & output$adj.P.Val < 0.05]
print(length(up.Genes))
print(length(dn.Genes))
GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE41621/GeneLst.RDS')

## DEG for met_express
#output <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE41621/Limma_output.txt',sep='\t',header = T)
sig.output <- output[abs(output$logFC)> log(1,2) & output$adj.P.Val<0.05,]
write.table(sig.output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE41621/DEG_P0.05_Lfc0.txt',sep='\t',row.names = F,quote = F)




##############################
### EMBO 2020
##############################
rm(list=ls())
library(readxl)
library(limma)

data <- read.csv('CRU/TALL/data/EMBO_2020/exprSet.txt',sep='\t',header  = T)
colnames(data)
label <- rep('tumor',ncol(data))
label[grep('^Thymu',colnames(data))]='normal'
Group <- data.frame(ID=colnames(data),ID_REF=label)

write.table(data,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/Gene_Expression_Value.txt',quote = F,sep='\t',row.names = T,col.names = T)
write.table(Group,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/DEG_design.txt',quote = F,sep='\t',row.names = F,col.names = F)


## 执行limma的DEG
print(max(data))

## 判断是否log转换
if(max(data) > 50){
  EXP=log(data+1,2)
}else{
  EXP=data
}
group_list <- factor(label,levels=c('normal','tumor'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) #FPKM

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene')  
print(output[output$gene %in% c('TAL1','S1PR3','ARRB1'),])
write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/Limma_output.txt',sep='\t',row.names = F,quote = F)

## DEGlst
up.Genes <- output$gene[output$logFC > 0 & output$adj.P.Val < 0.05]
dn.Genes <- output$gene[output$logFC < 0 & output$adj.P.Val < 0.05]
print(length(up.Genes))
print(length(dn.Genes))
GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/GeneLst.RDS')

## DEG for met_express
#output <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE41621/Limma_output.txt',sep='\t',header = T)
sig.output <- output[abs(output$logFC)> log(1,2) & output$adj.P.Val<0.05,]
write.table(sig.output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/DEG_P0.05_Lfc0.txt',sep='\t',row.names = F,quote = F)


##############################
### SY
##############################
rm(list=ls())
library(readxl)
library(limma)

data <- readRDS('/local/yanzijun/CRU/TALL/data/SY/RDS/RemoveBatchEffect_geneTPM.RDS')
colnames(data)

dele.samples <- c(
  # 'L434','L440','L445','L464','L488','L463','L449','L429','L484',
  #                 'L489','L494','L496','L495','L497','L498','L499','L500','L502',
  #                 'L503','L504','L505','L506','L507',
                  'L470','L459','L454','L245','L88')
dele.samples <- paste(dele.samples,'_TPM',sep='')
data <- data[,!colnames(data) %in% dele.samples]


Thy.samples <- colnames(data)[grep('^Thy',colnames(data))]
TALL.samples <- colnames(data)[grep('^L',colnames(data))]

data <- dplyr::select(data,c(Thy.samples,TALL.samples))
dim(data)
label <- c(rep('normal',length(Thy.samples)),rep('tumor',length(TALL.samples)))
Group <- data.frame(ID=colnames(data),ID_REF=label)


## 执行limma的DEG
print(max(data))

## 判断是否log转换
if(max(data) > 50){
  EXP=log(data+1,2)
}else{
  EXP=data
}
group_list <- factor(label,levels=c('normal','tumor'))
design <- model.matrix(~0+group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(EXP)

contrast.matrix <- makeContrasts(tumor-normal, levels=design)
fit <- lmFit(EXP, design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit,trend=T,robust=T) 

output <- topTable(fit, coef=1,n=Inf,adjust.method='fdr')
output <- tibble::rownames_to_column(output,'gene')  
print(output[output$gene %in% c('TAL1','S1PR3','ARRB1','PRMT5'),])

write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/Limma_output.txt',sep='\t',row.names = F,quote = F)

## DEGlst
up.Genes <- output$gene[output$logFC > 0 & output$adj.P.Val < 0.05]
dn.Genes <- output$gene[output$logFC < 0 & output$adj.P.Val < 0.05]
print(length(up.Genes))
print(length(dn.Genes))
GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/SY/GeneLst.RDS')

geneLst <- c('CDK6','TUBA1A','TUBB','TYMS','PSMB2')
print(output[output$gene %in% geneLst,])




df.other <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/Limma_output.txt',
                       header = TRUE,sep='\t')
print(df.other[df.other$gene %in% c('TAL1','S1PR3','ARRB1','PRMT5'),])

df.other <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/Limma_output.txt',
                       header = TRUE,sep='\t')
print(df.other[df.other$gene %in% c('TAL1','S1PR3','ARRB1','PRMT5'),])


df.other <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/Limma_output.txt',
                       header = TRUE,sep='\t')
print(df.other[df.other$gene %in% c('TAL1','S1PR3','ARRB1','PRMT5'),])



normal <- as.numeric(EXP[rownames(EXP)=='S1PR3',Thy.samples])
tumor <- as.numeric(EXP[rownames(EXP)=='S1PR3',TALL.samples])
S1PR3 <- data.frame(tissue=label,exp=c(normal,tumor))
wilcox.test(tumor,normal,alternative = 'less')

library(ggpubr)
library(ggsignif)
ggplot(S1PR3, aes(x=tissue, y=exp,fill=tissue)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1)+
  geom_jitter(width = 0.1) +
  theme_classic()+labs(title="S1PR3",y='The expression level',x='tissue')+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size=15,face="plain",color="black"),
        axis.title  = element_text(size=15,face="plain",color="black"),
        legend.position="none")+
  geom_signif(comparisons = list(c("normal", "tumor")),
              map_signif_level=TRUE)
