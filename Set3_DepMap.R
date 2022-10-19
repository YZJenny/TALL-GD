rm(list=ls())
library(ggplot2)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

### 提取出 Depmap 中 T-ALL细胞系中依赖性评分<0 的基因
#dep <- read.csv('/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_dependency.csv')
eff <- read.csv('/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect.csv')

META <- read.csv('/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/sample_info.csv')

eff <- eff[eff$DepMap_ID %in% intersect(eff$DepMap_ID,META$DepMap_ID),]
META <- META[META$DepMap_ID %in% intersect(eff$DepMap_ID,META$DepMap_ID),]
dim(eff);dim(META)

## T-ALL
TALL.celline <- META$DepMap_ID[META$Subtype=='Acute Lymphoblastic Leukemia (ALL), T-cell']
eff.TALL <- eff[eff$DepMap_ID %in% TALL.celline,]
rownames(eff.TALL) <- eff.TALL$DepMap_ID
eff.TALL$DepMap_ID <- NULL
dim(eff.TALL)
saveRDS(eff.TALL,'/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_TALL.RDS')

## B-ALL
BALL.celline <- META$DepMap_ID[META$Subtype=='Acute Lymphoblastic Leukemia (ALL), B-cell']
eff.BALL <- eff[eff$DepMap_ID %in% BALL.celline,]
rownames(eff.BALL) <- eff.BALL$DepMap_ID
eff.BALL$DepMap_ID <- NULL
dim(eff.BALL)
saveRDS(eff.BALL,'/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_BALL.RDS')


## AML
AML.celline <- META$DepMap_ID[grep('^Acute My',META$Subtype)]
eff.AML <- eff[eff$DepMap_ID %in% AML.celline,]
rownames(eff.AML) <- eff.AML$DepMap_ID
eff.AML$DepMap_ID <- NULL
dim(eff.AML)
saveRDS(eff.AML,'/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_AML.RDS')



## 所有细胞系都<0的基因
max.value <- apply(eff.TALL,2,max)

cutoff=0
if(cutoff==0){
  eff.Genes <- names(max.value)[max.value< cutoff]
}else{
  eff.Genes <- names(max.value)[max.value< -cutoff]
}
length(eff.Genes)

eff.Genes <- apply(as.matrix(eff.Genes),1,function(x) {unlist(strsplit(x,'\\.'))[1]})
saveRDS(eff.Genes,paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DepMap/effectGene_',cutoff,'.RDS',sep=''))

#### plot
eff.df <- eff.TALL[,colnames(eff.TALL) %in% eff.Genes]
dim(eff.df)
mean.value <- apply(eff.df,2,mean)
head(mean.value)
plot(density(mean.value))
boxplot(mean.value)
topN=50
top.genes <- names(sort(mean.value,decreasing = F)[1:topN])
top.df <- eff.df[,colnames(eff.df) %in% top.genes]
dim(top.df)
colnames(top.df) <- apply(as.matrix(top.genes),1,function(x) {unlist(strsplit(x,'\\.'))[1]})

p <- pheatmap(t(top.df),
              color = colorRampPalette(c("navy", "white"))(50),
              cluster_cols=F,cluster_rows = F, show_colnames = T,show_rownames = T,
              treeheight_col=0,
              fontsize = 15)
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/FIG/DepMap_Negative_top',topN,'.pdf',sep=''),p,width = 4.5,height = 10)


### barplot
get_freq <- function(x){
  freq=length(which(x<0))
  return(freq)
}

Neg.Freq <- as.data.frame(apply(eff.TALL,2,function(x) get_freq(x)))
colnames(Neg.Freq)[1] <- "Class" 
Neg.Freq$Class <- factor(Neg.Freq$Class,levels = c(5,4,3,2,1,0))
Neg.Freq <- as.data.frame(table(Neg.Freq$Class))
Neg.Freq <- data.frame(Class=Neg.Freq$Var1,Freq=c(cumsum(Neg.Freq$Freq[1:5]),Neg.Freq$Freq[6]))

Neg.Freq$Group='=5'
Neg.Freq$Group[Neg.Freq$Class=='4']='>=4'
Neg.Freq$Group[Neg.Freq$Class=='3']='>=3'
Neg.Freq$Group[Neg.Freq$Class=='2']='>=2'
Neg.Freq$Group[Neg.Freq$Class=='1']='>=1'
Neg.Freq$Group[Neg.Freq$Class=='0']='=0'
Neg.Freq$Group=factor(Neg.Freq$Group,levels=c('=5','>=4','>=3','>=2','>=1','=0'))

Neg.Freq$type='B'
Neg.Freq$type[Neg.Freq$Class=='5']='A'


p <-ggplot(data =Neg.Freq, aes(x = Group, y = Freq,
                            fill = type)) + geom_bar(stat = "identity")+
  theme_classic2()+
  geom_text(aes(label = Freq), vjust = -0.8)+
  scale_fill_manual(values = c(getPalette[1],"#999999"))+
  scale_y_continuous(expand = c(0,0))+ # 设置y轴从0开始
  theme(panel.background = element_blank(), # 去掉背景格子
        axis.line.x = element_line(colour = "black"),
        axis.text = element_text(colour = 'black'),
        legend.position = 'none')
ggsave('/local/yanzijun/CRU/TALL_FM/res/FIG/DepMap_Freq.pdf',p,width = 4,height = 4)
