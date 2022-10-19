#### Combat batch effect
rm(list=ls())
library(mgcv)
library(nlme)
library(genefilter)
library(BiocParallel)
library(sva)
library(dplyr)
library(edgeR)

zero.rows.del <- function(data,batch){
  data <- as.matrix(data)
  batch <- as.factor(batch)
  zero.rows.lst <- lapply(levels(batch), function(batch_level) {
    if (sum(batch == batch_level) > 1) {
      return(which(apply(data[, batch == batch_level], 1,
                         function(x) {
                           var(x) == 0
                         })))
    }else {
      
      return(which(rep(1, 3) == 2))
      
    }       
  })
  zero.rows <- Reduce(union, zero.rows.lst)
  keep.rows <- setdiff(1:nrow(data), zero.rows)
  if (length(zero.rows) > 0) {
    cat(sprintf("Found %d genes with uniform expression within a single batch (all zeros); these will not be adjusted for batch.\n",
                length(zero.rows)))
    data.orig <- data
    data <- data[keep.rows, ]
  }
  return(data)
}

EXP0=read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/Gene_Expression_Value.txt',header=T,sep='\t')
EXP0 <- tibble::rownames_to_column(EXP0,'Gene')
Group0=read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/EMBO/DEG_design.txt',header=F,sep='\t')
Label0=Group0[,2]

EXP1=read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/Gene_Expression_Value.txt',header=T,sep='\t')
EXP1 <- tibble::rownames_to_column(EXP1,'Gene')
Group1=read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE26713/DEG_design.txt',header=F,sep='\t')
Label1=Group1[,2]

EXP2=read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/Gene_Expression_Value.txt',header=T,sep='\t')
EXP2 <- tibble::rownames_to_column(EXP2,'Gene')
Group2=read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/GSE48558/DEG_design.txt',header=F,sep='\t')
Label2=Group2[,2]


data <- merge(merge(EXP0,EXP1,by='Gene'),EXP2,by='Gene')
rownames(data) <- data$Gene
data$Gene <- NULL
dim(data) 

label <- factor(c(Label0,Label1,Label2),levels = c('tumor','normal'))

batch <- c(rep('EMBO',length(Label0)),rep('GSE26713',length(Label1)),rep('GSE48558',length(Label2)))

data.no0 <- zero.rows.del(data,batch)

## 执行limma的removeBatchEffect
# library(limma)
# expr_batch <- removeBatchEffect(data.no0, batch)

## 执行Combat
mod <- model.matrix(~label)
expr_batch <- ComBat(dat = data.no0, batch = batch, mod = mod,mean.only=TRUE)
dim(expr_batch)

combat.res <- list(EXP=expr_batch, Group=label, Batch=batch)
saveRDS(combat.res,'/local/yanzijun/CRU/TALL_FM/res/Combat/combat_res.RDS')

########################
## 1. 检验效果
########################
library("FactoMineR")
library("ggplot2")
library("factoextra")
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

pca.plot = function(dat,col){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               #addEllipses = TRUE,
               legend.title = "Group"
  )
}


p1 <- pca.plot(data,factor(batch))+scale_color_manual(values = getPalette(9)[1:3])## 处理前
p2 <- pca.plot(expr_batch,factor(batch))+scale_color_manual(values = getPalette(9)[1:3])  ## 处理后

p3 <- pca.plot(data,factor(label))+scale_color_manual(values = getPalette(9)[1:2])## 处理前
p4 <- pca.plot(expr_batch,factor(label))+scale_color_manual(values = getPalette(9)[1:2])  ## 处理后

PCA <- ggarrange(p1,p2,p3,p4,nrow = 2,ncol = 2)
pdf('/local/yanzijun/CRU/TALL_FM/res/FIG/Combat_PCA.pdf',width = 8,height = 6)
print(PCA)
dev.off()



########################
## 2. PCA
########################
rm(list=ls())
.libPaths(c( "/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
             "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library("FactoMineR")
library("ggplot2")
library("factoextra")
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

combat.res <- readRDS('/local/yanzijun/CRU/TALL_FM/res/Combat/combat_res.RDS')
EXP <- combat.res$EXP
Group <- combat.res$Group
exp <- EXP[rownames(EXP)%in%candiGenes,]
PCA <- pca.plot(exp,factor(Group))+
  scale_color_manual(values = getPalette(9)[1:2])
print(PCA)

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/PCA_',cutoff,'_Combat.pdf',sep=''),width = 4,height = 3)
print(PCA)
dev.off()


########################
## 3. correlation 相关性热图
########################
rm(list=ls())
library(ggplot2)
library(ggsci)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)


# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  #cormat[lower.tri(cormat)]<- NA
  cormat[upper.tri(cormat)]<- NA
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
  #cormat <- round(cor(t(data)),2)
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = getPalette[2], high = getPalette[1], mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1,colour = 'black'),
          axis.text.y = element_text(size = 9, colour = 'black'))+
    coord_fixed()
  p <- ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 1) +
    theme(axis.title = element_blank(),axis.text  = element_text(colour = 'black'),
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

combat.res <- readRDS('/local/yanzijun/CRU/TALL_FM/res/Combat/combat_res.RDS')
EXP <- combat.res$EXP
Group <- combat.res$Group


case=EXP[rownames(EXP)%in%candiGenes,which(Group=='tumor')]
ctrl=EXP[rownames(EXP)%in%candiGenes,which(Group=='normal')]

p.tumor <- plot_cor(data=case)
p.normal <- plot_cor(data=ctrl)
p.cor <- ggarrange(p.tumor,p.normal,nrow = 1)
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/Cor_',cutoff,'_Combat_up.pdf',sep=''),p.cor,width = 12,height = 5)
