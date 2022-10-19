################## 
#### 计算ATAC vs RNA vs DeapMap cor
##################
rm(list=ls())
ATAC_output <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/Limma_EV7_gene.RDS')
## 选取logFC最大
ATAC_output <- ATAC_output[order(ATAC_output$gene,ATAC_output$logFC,decreasing=T),]
ATAC_output=ATAC_output[!duplicated(ATAC_output$gene),]
colnames(ATAC_output)[2] <- 'logFC_ATAC'


GSE_ID='EMBO'
RNA_output <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/',GSE_ID,'/Limma_output.txt',sep=''),header = T,sep='\t')
RNA_output <- RNA_output[,c('gene','logFC')]
colnames(RNA_output)[2] <- 'logFC_RNA'

### 整体基因的cor
plot.df <- merge(RNA_output,ATAC_output,by='gene')
cor.test(plot.df$logFC_RNA,plot.df$logFC_ATAC)

### 只展示29候选基因的 cor
candiGenes <- readRDS('/local/yanzijun/CRU/TALL_FM/res/candiSet4/candiGenes_0.RDS')
length(candiGenes)
subplot.df <- plot.df[plot.df$gene %in% candiGenes,]
dim(subplot.df)
res <- cor.test(subplot.df$logFC_RNA,subplot.df$logFC_ATAC)
res$p.value
res$estimate

library(ggrepel)
p.cor <- ggplot(data=subplot.df, aes(x=logFC_RNA, y=logFC_ATAC)) + 
  geom_point(alpha=0.7,color = 'red') +stat_smooth(method="lm")+
  geom_text_repel(aes(logFC_RNA, logFC_ATAC, label = gene))+
  theme_classic(base_size = 10)+
  geom_text(aes(0.1, 1.2, label = paste("Cor=",round(res$estimate,2),'\n',
                                        "p-value=",round(res$p.value,3),sep='')))

ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/candiSet4/FIG/CorPlot_',GSE_ID,'.pdf',sep=''),p.cor,
       width = 6,height = 5  )


### 三个数据集取均值
RNA_output1 <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/Limma_output.txt',header = T,sep='\t')
RNA_output1 <- RNA_output1[,c('gene','logFC')]
RNA_output2 <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/Limma_output.txt',header = T,sep='\t')
RNA_output2 <- RNA_output2[,c('gene','logFC')]
RNA_output3 <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/Limma_output.txt',header = T,sep='\t')
RNA_output3 <- RNA_output3[,c('gene','logFC')]

RNA_output <- merge(merge(RNA_output1,RNA_output2,by='gene'),RNA_output3,'gene')
dim(RNA_output)
RNA_output$logFC_RNA <- apply(RNA_output[,-1],1,mean)
RNA_output <- RNA_output[,c('gene','logFC_RNA')]
plot.df <- merge(RNA_output,ATAC_output,by='gene')
dim(plot.df)
cor.test(plot.df$logFC_RNA,plot.df$logFC_ATAC) # 0.1

subplot.df <- plot.df[plot.df$gene %in% candiGenes,]
dim(subplot.df)
res <- cor.test(subplot.df$logFC_RNA,subplot.df$logFC_ATAC)
res$p.value
res$estimate # 0.26

########################
## combat后的数据做limma
########################
library(limma)
rm(list=ls())
plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

combat.res <- readRDS('/local/yanzijun/CRU/TALL_FM/res/Combat/combat_res.RDS')
data <- combat.res$EXP
label <- combat.res$Group


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
dim(output)
output[1:2,]
saveRDS(output,'/local/yanzijun/CRU/TALL_FM/res/Combat/limma_res.RDS')

################## 
#### 计算ATAC vs RNA vs DeapMap cor
##################
ATAC_output <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/Limma_EV7_gene.RDS')
## 选取logFC最大
ATAC_output <- ATAC_output[order(ATAC_output$gene,ATAC_output$logFC,decreasing=T),]
ATAC_output=ATAC_output[!duplicated(ATAC_output$gene),]
colnames(ATAC_output)[2] <- 'logFC_ATAC'

RNA_output <- readRDS('/local/yanzijun/CRU/TALL_FM/res/Combat/limma_res.RDS')
RNA_output <- RNA_output[,c('gene','logFC')]
colnames(RNA_output)[2] <- 'logFC_RNA'

DepMap <- readRDS('/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_TALL.RDS')
colnames(DepMap) <- apply(as.matrix(colnames(DepMap)),1,function(x) {unlist(strsplit(x,'\\.'))[1]})
DepMap <- data.frame(gene=colnames(DepMap),DepMap=apply(DepMap,2,mean))
DepMap[1:2,]

candiGenes <- readRDS('/local/yanzijun/CRU/TALL_FM/res/candiSet4/candiGenes_0.RDS')
length(candiGenes)


### 1. RNA vs ATAC 的cor
plot.df <- merge(RNA_output,ATAC_output,by='gene') #整体基因
cor.test(plot.df$logFC_RNA,plot.df$logFC_ATAC) #0.07

### 只展示29候选基因的 cor
subplot.df <- plot.df[plot.df$gene %in% candiGenes,]
dim(subplot.df)
res <- cor.test(subplot.df$logFC_RNA,subplot.df$logFC_ATAC,method = 'spearman')
res$p.value
res$estimate

library(ggrepel)
p.cor <- ggplot(data=subplot.df, aes(x=logFC_RNA, y=logFC_ATAC)) + 
  geom_point(alpha=0.7,color = 'red') +stat_smooth(method="lm")+
  geom_text_repel(aes(logFC_RNA, logFC_ATAC, label = gene))+
  theme_classic(base_size = 10)+
  geom_text(aes(0.5, 1.2, label = paste("Cor=",round(res$estimate,2),'\n',
                                        "p-value=",round(res$p.value,3),sep='')))

ggsave('/local/yanzijun/CRU/TALL_FM/res/candiSet4/FIG/CorPlot_combat_RNA_ATAC.pdf',p.cor,
       width = 6,height = 5)


### 2. RNA vs DepMap 的cor
plot.df <- merge(RNA_output,DepMap,by='gene') #整体基因
cor.test(plot.df$logFC_RNA,plot.df$DepMap,method = 'spearman') # -0.14

### 只展示29候选基因的 cor
subplot.df <- plot.df[plot.df$gene %in% candiGenes,]
dim(subplot.df)
res <- cor.test(subplot.df$logFC_RNA,subplot.df$DepMap,method = 'spearman')
res$p.value
res$estimate


### 3. ATAC vs DepMap 的cor
plot.df <- merge(ATAC_output,DepMap,by='gene') #整体基因
cor.test(plot.df$logFC_ATAC,plot.df$DepMap,method = 'spearman') # -0.14

### 只展示29候选基因的 cor
subplot.df <- plot.df[plot.df$gene %in% candiGenes,]
dim(subplot.df)
res <- cor.test(subplot.df$logFC_ATAC,subplot.df$DepMap,method = 'spearman')
res$p.value
res$estimate
