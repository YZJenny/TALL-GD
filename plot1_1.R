rm(list=ls())
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)
library(pheatmap)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)


############  
### 1. volcaono, heatmap and vlnplot/pheatmap
############
## 1.1 volcaono
file.list <- c('GSE26713','GSE48558','EMBO')
for(file in file.list){
  print(file)
  
  limma_res <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/',file,'/Limma_output.txt',sep=''),header = T,sep='\t') 
  data <- limma_res
  
  qvalue <- 0.05
  fc_cutoff <- 1
  data <- mutate(data, 
                 sig=ifelse((data$adj.P.Val < qvalue & data$logFC > log(fc_cutoff,2))| 
                              (data$adj.P.Val < qvalue & data$logFC < -log(fc_cutoff,2)) ,
                            ifelse(data$logFC > log(fc_cutoff,2),'UP','DOWN'),'no'))
  data <- na.omit(data)
  data$log10Qvalue=-log10(data$adj.P.Val)
  
  figure <- ggplot(data,aes(x=logFC,y=log10Qvalue))+geom_point(aes(color=sig))+
    scale_color_manual(values = c(mycol[2],"#999999",mycol[1]))+
    theme_bw()+
    labs(title=file,x="log2 (Fold Change)",y = "-log10 (adjP value)")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=12,colour = 'black'),
          axis.line = element_line(size=0.5),
          panel.grid.major=element_line(colour=NA),panel.grid=element_blank())
  p <- figure+
    geom_vline(xintercept=log(fc_cutoff,2),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(qvalue),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
  ggsave(plot=p,filename = paste('/local/yanzijun/CRU/TALL_FM/res/FIG/volcano_',file,'.pdf',sep=''),
         width = 4,height = 3.5)
}




### 1.2 heatmap
file.list <- c('GSE26713','GSE48558','EMBO')
for(file in file.list){
  print(file)
  
  limma_res <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/',file,'/Limma_output.txt',sep=''),header = T,sep='\t') 
  EXP=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/Gene_Expression_Value.txt',sep=''),header=T,sep='\t')
  Group=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/DEG_design.txt',sep=''),header=F,sep='\t')
  Label=Group[,2]

  up.sigG <- limma_res[which(limma_res$adj.P.Val < 0.05 & limma_res$logFC > log(1,2) ),]
  up.sigG <- up.sigG$gene[order(up.sigG$logFC,decreasing = T)][1:30]

  dn.sigG <- limma_res[which(limma_res$adj.P.Val < 0.05 & limma_res$logFC < -log(1,2) ),]
  dn.sigG <- dn.sigG$gene[order(dn.sigG$logFC,decreasing =F )][1:30]
  sigG <- c(up.sigG,dn.sigG)
  
  #sigG <- candiGenes
  print(length(sigG))
  
  ## heatmap
  df_heatmap <- EXP[which(rownames(EXP)%in% sigG),]
  print(dim(df_heatmap))
  
  annotation_col = data.frame(SampleType = factor(Label))
  rownames(annotation_col) = colnames(df_heatmap)
  
  p <- pheatmap(df_heatmap, annotation_col = annotation_col,scale = 'row',
                color = colorRampPalette(c("navy", "white", "firebrick2"))(50),
                cluster_cols=T,cluster_rows = T, show_colnames = F,show_rownames = F,
                treeheight_col=0,
                breaks = seq(-2, 2, length.out = 50),fontsize = 15)
  ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/FIG/heatmap_',file,'.pdf',sep=''),p,width = 7,height = 6)
}


## 1.3 vlnplot
plan='candiSet4'
cutoff='0' 

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

file.list <- c('GSE26713','GSE48558','EMBO')
for(file in file.list){
  print(file)
  
  EXP=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/Gene_Expression_Value.txt',sep=''),header=T,sep='\t')
  Group=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/DEG_design.txt',sep=''),header=F,sep='\t')
  Label=Group[,2]
  print(table(Label))

  topN=length(candiGenes)
  sub.EXP <- EXP[rownames(EXP) %in% candiGenes[1:topN],]
  plot.df <- reshape2::melt(t(sub.EXP))
  head(plot.df)
  plot.df <- merge(plot.df,Group,by.x='Var1',by.y = 'V1')
  colnames(plot.df) <- c('sample','gene','exp','group')
  plot.df$group=factor(plot.df$group,levels = c('tumor','normal'))
  #plot.df$gene <- factor(plot.df$gene,levels = candiGenes[1:topN])
  
  p.vln <- ggboxplot(plot.df, x = "group", y = "exp",
                     color = "group", 
                     add = NULL,
                     facet.by = "gene", short.panel.labs = FALSE,nrow =3)+
    stat_compare_means(label = "p.format")+
    scale_color_manual(values = mycol[1:2])+
    scale_fill_manual(values = mycol[1:2])
  ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/vlnplot_',file,'_',cutoff,'.pdf',sep=''),
         p.vln,width = 12,height = 7)
}


################
#### pheatmap 
################
file.list <- c('GSE26713','GSE48558','EMBO')
for(file in file.list){
  print(file)
  
  EXP=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/Gene_Expression_Value.txt',sep=''),header=T,sep='\t')
  Group=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/DEG_design.txt',sep=''),header=F,sep='\t')
  rownames(Group)=Group$V1
  
  normal=Group$V1[Group$V2=='normal']
  tumor=Group$V1[Group$V2=='tumor']
  Group <- as.data.frame(t(select(as.data.frame(t(Group)),c(normal,tumor))))
  Group$V1=NULL;colnames(Group)[1]='Group'
  head(Group)
  
  sub.EXP <- EXP[rownames(EXP) %in% candiGenes,]
  sub.EXP <- select(sub.EXP,rownames(Group))
  print(all(colnames(sub.EXP)==rownames(Group)))
  exp.plot <- scale(sub.EXP)
  
  ann_colors = list(
    Group=c(normal=mycol[2],tumor=mycol[1])
  )
  
  p <- pheatmap::pheatmap(exp.plot,cluster_cols = F,cluster_rows = T, 
                          annotation_col = Group,show_colnames = F,
                          treeheight_row=0,treeheight_col=0,
                          annotation_colors = ann_colors,border_color = "white",
                          colorRampPalette(c("navy","white", "firebrick3"))(50),
                          gaps_col = cumsum(c(length(normal),length(tumor)))
  )
  p
  
  ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/Heatmap_',file,'_',cutoff,'.pdf',sep=''),
         p,width = 10,height = 7)
}



##############################
### 2. dependency score 在 白血病(T-ALL/B-ALL/AML)的cell line中特异性比较
##############################
rm(list=ls())
library(ggplot2)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))

get_plotdf <- function(candiGenes,eff.input,group,cutoff){
  eff.df <- readRDS(eff.input)
  dim(eff.df)
  
  col.Genes <- apply(as.matrix(colnames(eff.df)),1,function(x) {unlist(strsplit(x,'\\.'))[1]})
  colnames(eff.df) <-col.Genes 
  df <- eff.df[,colnames(eff.df)%in%candiGenes]
  dim(df)
  med <- apply(df,2,median)
  rank <- names(sort(med,decreasing = F))
  
  plot.df <- reshape2::melt(df)
  
  if(cutoff=='0'|cutoff=='0.5'){
    plot.df$value <- abs(plot.df$value)
  }
  colnames(plot.df) <- c('gene','abs.EffectScore')
  plot.df$group <- group
  plot.df$gene=factor(plot.df$gene,levels = rank)
  return(plot.df)
}



plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

## 在 TALL的cell line中
plot.TALL <- get_plotdf(candiGenes=candiGenes,
                        eff.input='/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_TALL.RDS',
                        group = 'TALL',
                        cutoff=cutoff)
p.TALL <- ggplot(data = plot.TALL) + 
  geom_boxplot(aes(x = gene, y = abs.EffectScore, color = factor(gene))) + 
  geom_jitter(aes(x = gene, y = abs.EffectScore, color = factor(gene)), 
              position = position_jitterdodge())+
  theme_classic()+
  theme(legend.position ='none',
        axis.text.x = element_text(hjust = 1,angle = 30))+
  scale_color_manual(values = mycol(length(candiGenes)))

if(cutoff=='surv'){
  p.TALL <- p.TALL+
    geom_hline(aes(yintercept=0),colour="black", linetype="dashed")
}

p.TALL
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/EffectScore_TALL_',cutoff,'.pdf',sep=''),width = 7,height = 4)
print(p.TALL)
dev.off()

### flip plot
plot.TALL.flip <- plot.TALL
plot.TALL.flip$gene <- factor(plot.TALL.flip$gene,levels = rev(levels(plot.TALL$gene)))

p.TALL.flip <- ggplot(data = plot.TALL.flip) + 
  geom_boxplot(aes(x = gene, y = abs.EffectScore,
                   fill = factor(gene)
                   )) + 
  geom_jitter(aes(x = gene, y = abs.EffectScore,
                  fill = factor(gene)), 
              position = position_jitterdodge())+
  scale_fill_manual(values = mycol(length(candiGenes)))+
  scale_color_manual(values = mycol(length(candiGenes)))+
  theme_classic()+
  theme(legend.position ='none',
        axis.text = element_text(color = 'black')
        #axis.text.x = element_text(hjust = 1,angle = 30)
        )+
  coord_flip()
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/EffectScore_TALL_',cutoff,'_flip.pdf',sep=''),width = 4,height = 7)
print(p.TALL.flip)
dev.off()




## 在 白血病(T-ALL/B-ALL/AML)
plot.BALL <- get_plotdf(candiGenes = candiGenes,
                        eff.input = '/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_BALL.RDS',
                        group = 'BALL',
                        cutoff = cutoff)
plot.AML <- get_plotdf(candiGenes=candiGenes,
                       eff.input='/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CRISPR_gene_effect_AML.RDS',
                       group='AML',
                       cutoff = cutoff)

plot.L <- rbind(plot.TALL,plot.BALL,plot.AML)
plot.L$group=factor(plot.L$group,levels = c('TALL','BALL','AML'))
# p <- ggplot(data = plot.df) + 
#   geom_boxplot(aes(x = gene, y = abs.EffectScore, color = group)) + 
#   geom_jitter(aes(x = gene, y = abs.EffectScore, color = group), 
#               position = position_jitterdodge())+
#   theme_classic()

p.L <- ggboxplot(plot.L, x = "group", y = "abs.EffectScore",
               color = "group", palette = "jco",
               add = "jitter",
               facet.by = "gene", short.panel.labs = FALSE,nrow =2)+
  stat_compare_means(label = "p.format")


plot.ALL <- rbind(plot.TALL,plot.BALL)
plot.ALL$group=factor(plot.ALL$group,levels = c('TALL','BALL'))
p.ALL <- ggboxplot(plot.ALL, x = "group", y = "abs.EffectScore",
                 color = "group", palette = "jco",
                 add = "jitter",
                 facet.by = "gene", short.panel.labs = FALSE,nrow =2)+
  stat_compare_means(label = "p.format")


plot.LM <- rbind(plot.TALL,plot.AML)
plot.LM$group=factor(plot.LM$group,levels = c('TALL','AML'))
p.LM <- ggboxplot(plot.LM, x = "group", y = "abs.EffectScore",
                   color = "group", palette = "jco",
                   add = "jitter",
                   facet.by = "gene", short.panel.labs = FALSE,nrow =2)+
  stat_compare_means(label = "p.format")

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/EffectScore_',cutoff,'.pdf',sep=''),width = 16,height = 5)
print(p.LM)
print(p.ALL)
print(p.L)
dev.off()
### 结论：effectscore在TALL/BALL/AML中没有差异



##############################
### 2. gene expression的PCA分析
##############################
rm(list=ls())
.libPaths(c( "/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
             "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library("FactoMineR")
library("ggplot2")
library("factoextra")
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))

pca.plot = function(dat,col){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               #addEllipses = TRUE,
               legend.title = "Groups"
  )
}


plan='candiSet2'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

## EMBO
EXP <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/Gene_Expression_Value.txt',sep='\t',header = T )
Group <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/DEG_design.txt',sep='\t',header = F)
Group <- factor(Group$V2,levels = c('tumor','normal'))
exp <- EXP[rownames(EXP)%in%candiGenes,]
PCA.EMBO <- pca.plot(exp,factor(Group))+
  scale_color_manual(values = mycol(9)[1:2])
print(PCA.EMBO)

## GSE26713
EXP <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/Gene_Expression_Value.txt',sep='\t',header = T )
Group <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/DEG_design.txt',sep='\t',header = F)
Group <- factor(Group$V2,levels = c('tumor','normal'))
exp <- EXP[rownames(EXP)%in%candiGenes,]
PCA.GSE1 <- pca.plot(exp,col=factor(Group))+
  scale_color_manual(values = mycol(9)[1:2])
print(PCA.GSE1)

## GSE48558
EXP <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/Gene_Expression_Value.txt',sep='\t',header = T )
Group <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/DEG_design.txt',sep='\t',header = F)
Group <- factor(Group$V2,levels = c('tumor','normal'))
exp <- EXP[rownames(EXP)%in%candiGenes,]
PCA.GSE2 <- pca.plot(exp,factor(Group))+
  scale_color_manual(values = mycol(9)[1:2])
print(PCA.GSE2)


### TARGET T-ALL VS AML
library(readxl)
library(dplyr)

## AML
AML_exp <- read.table('/local/yanzijun/public/TARGET/AML/TARGET_NBM_AML_QuantileNormalized_RPKM.txt',header = T)
AML_exp <- AML_exp[!duplicated(AML_exp$gene_name),]
rownames(AML_exp) <- AML_exp$gene_name;AML_exp$gene_name=AML_exp$gene_id=NULL

AML_exp <- AML_exp[,-grep('^TARGET.00',colnames(AML_exp))] #删除20个正常BM样本
AML <- AML_exp[rownames(AML_exp)%in%candiGenes,]
AML <- AML[order(rownames(AML)),]

## B-ALL
BALL_exp <- read.csv('/local/yanzijun/public/cBioportal/all_phase2_target_2018_pub/data_mrna_seq_rpkm.txt',sep='\t',header = T)
BALL_exp <- BALL_exp[!duplicated(BALL_exp$Hugo_Symbol),]
rownames(BALL_exp) <- BALL_exp$Hugo_Symbol;BALL_exp$Hugo_Symbol=BALL_exp$Entrez_Gene_Id=NULL

BALL_clin <- read.csv('/local/yanzijun/public/cBioportal/all_phase2_target_2018_pub/data_clinical_sample.txt',sep='\t',skip = 4,header = T)
table(BALL_clin$CELL_OF_ORIGIN)
table(BALL_clin$CANCER_TYPE_DETAILED)

BALL <- BALL_exp[rownames(BALL_exp)%in% candiGenes,]
BALL <- BALL[order(rownames(BALL)),]
dim(BALL)

## ALAL
ALAL_exp <- read.csv('/local/yanzijun/public/TARGET/ALAL/MPAL-allGeneExpr-90S-rlog.txt',header = T,sep='\t')
ALAL_exp <- ALAL_exp[!duplicated(ALAL_exp$Gene),]
rownames(ALAL_exp) <- ALAL_exp$Gene;ALAL_exp$Gene=ALAL_exp$geneIDs=NULL

ALAL <- ALAL_exp[rownames(ALAL_exp)%in%candiGenes,]
ALAL <- ALAL[order(rownames(ALAL)),]


## T-ALL
library(readxl)
TALL_exp <- as.data.frame(read_excel("/local/yanzijun/CRU/TALL/data/NG_2017/RNAseq_FPKM.xlsx"))
rownames(TALL_exp) <- TALL_exp$Gene;TALL_exp$Gene <- NULL
TALL <- TALL_exp[rownames(TALL_exp)%in%candiGenes,]
TALL <- TALL[order(rownames(TALL)),]
dim(TALL)

## merge
print(all(rownames(TALL)==rownames(AML)))
print(all(rownames(TALL)==rownames(BALL)))
print(all(rownames(TALL)==rownames(ALAL)))

ol.genes <- intersect(rownames(AML),intersect(rownames(TALL),intersect(rownames(BALL),rownames(ALAL))))
length(ol.genes)
TALL <- TALL[rownames(TALL)%in%ol.genes,]
TALL <- TALL[order(rownames(TALL)),]

AML <- AML[rownames(AML)%in%ol.genes,]
AML <- AML[order(rownames(AML)),]

BALL <- BALL[rownames(BALL)%in%ol.genes,]
BALL <- BALL[order(rownames(BALL)),]

ALAL <- ALAL[rownames(ALAL)%in%ol.genes,]
ALAL <- ALAL[order(rownames(ALAL)),]

print(all(rownames(TALL)==rownames(AML)))
print(all(rownames(TALL)==rownames(BALL)))
print(all(rownames(TALL)==rownames(ALAL)))


TALL <- tibble::rownames_to_column(TALL,'gene')
AML <- tibble::rownames_to_column(AML,'gene')
BALL <- tibble::rownames_to_column(BALL,'gene')
ALAL <- tibble::rownames_to_column(ALAL,'gene')
EXP <- merge(merge(merge(TALL,AML,by='gene'),BALL,by='gene'),ALAL,by='gene')
EXP[1:3,1:3]
rownames(EXP) <- EXP$gene;EXP$gene=NULL
#EXP <- cbind(cbind(cbind(TALL,AML),BALL),ALAL)
Group <- c(rep('T-ALL',ncol(TALL)-1),rep('AML',ncol(AML)-1),rep('B-ALL',ncol(BALL)-1),rep('ALAL',ncol(ALAL)-1))
Group <- factor(Group,levels = c('T-ALL','AML','B-ALL','ALAL'))


PCA.L <- pca.plot(EXP,factor(Group))+
  scale_color_manual(values = mycol(9)[1:4])
print(PCA.L)


PCA <- ggarrange(PCA.GSE1,PCA.GSE2,PCA.EMBO,PCA.L,nrow = 2,ncol = 2)
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/PCA_',cutoff,'.pdf',sep=''),width = 8,height = 6)
print(PCA)
dev.off()


#####################
### 4. correlation图
#####################
rm(list=ls())
library(ggplot2)
library(ggsci)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)


# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  return(cormat)
}

plot_cor <- function(data){
  cormat <- round(cor(t(data)),2)
  #cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = mycol[2], high = mycol[1], mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1))+
    coord_fixed()
  p <- ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 1) +
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.grid.major = element_blank(),panel.border = element_blank(),
          panel.background = element_blank(),axis.ticks = element_blank(),
          legend.justification = c(1, 0),legend.position = c(0.6, 0.7),legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
  return(p)
}


plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

file.list=c('GSE26713','GSE48558','EMBO')

for(file in file.list){
  print(file)
  EXP=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/Gene_Expression_Value.txt',sep=''),header=T,sep='\t')
  Group=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/DEG_design.txt',sep=''),header=F,sep='\t')
  
  case=EXP[rownames(EXP)%in%candiGenes,Group$V1[Group$V2=='tumor']]
  ctrl=EXP[rownames(EXP)%in%candiGenes,Group$V1[Group$V2=='normal']]
  
  p.tumor <- plot_cor(data=case)
  p.normal <- plot_cor(data=ctrl)
  p.cor <- ggarrange(p.tumor,p.normal,nrow = 1)
  ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/Cor_',cutoff,'_',file,'.pdf',sep=''),p.cor,width = 12,height = 5)
}


#####################
## 不同数据集DEG的FC的correlation图, cor值都比较小
#####################
rm(list=ls())
library(ggplot2)

DEG0 <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/Limma_output.txt',header = T)
DEG0.gene <- DEG0$gene[DEG0$logFC>0 & DEG0$adj.P.Val<0.05]
DEG1 <- read.csv('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/Limma_output.txt',header = T,sep='\t')
DEG1.gene <- DEG1$gene[DEG1$logFC>0 & DEG1$adj.P.Val<0.05]
DEG2 <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/Limma_output.txt',header = T)
DEG2.gene <- DEG2$gene[DEG2$logFC>0 & DEG2$adj.P.Val<0.05]

ol.genes <- intersect(intersect(DEG0.gene,DEG1.gene),DEG2.gene)
print(length(ol.genes))

data <- merge(DEG1,DEG2,by='gene')
data <- data[data$gene %in% ol.genes,c('gene','logFC.x','logFC.y')]
dim(data)
cor.test(data$logFC.x,data$logFC.y,method = 'spearman')
ggplot(data = data,aes(x=logFC.x,y=logFC.y))+
  geom_point()
