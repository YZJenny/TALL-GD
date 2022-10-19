#### 29个基因的heatmap
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

plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)


### heatmap
file.list <- c('GSE26713','GSE48558','EMBO','Combat')
for(file in file.list){
  print(file)
  
  if(file!='Combat'){
    limma_res <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/',file,'/Limma_output.txt',sep=''),header = T,sep='\t') 
    EXP=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/Gene_Expression_Value.txt',sep=''),header=T,sep='\t')
    Group=read.table(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/met_express/',file,'/DEG_design.txt',sep=''),header=F,sep='\t')
    Label=Group[,2]
  }else{
    limma_res <- readRDS('/local/yanzijun/CRU/TALL_FM/res/Combat/limma_res.RDS')
    combat.res <- readRDS('/local/yanzijun/CRU/TALL_FM/res/Combat/combat_res.RDS')
    EXP <- as.data.frame(combat.res$EXP)
    Label <- combat.res$Group
  }

  sample_order <- c(colnames(EXP)[which(Label=='normal')],
                    colnames(EXP)[which(Label=='tumor')])
  EXP <- select(EXP,sample_order)
  
  sigG <- candiGenes
  print(length(sigG))
  
  ## heatmap
  df_heatmap <- EXP[which(rownames(EXP)%in% sigG),]
  print(dim(df_heatmap))
  
  annotation_col = data.frame(SampleType = c(rep('normal',length(which(Label=='normal'))),
                                             rep('tumor',length(which(Label=='tumor')))))
  rownames(annotation_col) = colnames(df_heatmap)
  
  ann_colors = list(
    SampleType=c(normal=mycol[2],tumor=mycol[1]))
  
  
  p.heatmap <- pheatmap::pheatmap(df_heatmap,scale = 'row',cluster_cols = T,cluster_rows = T, 
                                  clustering_method = 'median',
                                  annotation_col = annotation_col,show_colnames = F,show_rownames = T,
                                  treeheight_row=0,treeheight_col=12,
                                  annotation_colors = ann_colors,border_color = "white",
                                  colorRampPalette(c("navy","white",  "firebrick3"))(length(seq(-2,2,by = 0.1))),
                                  gaps_col = cumsum(
                                    c(length(which(Label=='normal')),length(which(Label=='tumor'))) ),
                                  breaks = seq(-2,2,by = 0.1),legend_breaks = seq(-2,2,1))
  ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/candiSet4/FIG/heatmap_',file,'_candiGene.pdf',sep=''),
         p.heatmap,width = 10,height = 5)
}
