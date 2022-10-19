rm(list=ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

plan='candiSet4'
cutoff='0'

setwd(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/PPI',sep=''))

KEGG <- read.table('enrichment.KEGG.tsv',header = F,sep='\t')
plot.df <- KEGG
colnames(plot.df) <- c('pathway','Description','Counts','AllCounts','Strength','FDR','ENSG','Symbol')
print(plot.df$Description)
#plot.df <- plot.df[c(2,3,4,7,11,13),]
plot.df$log10Qval <- -log(plot.df$FDR,10)
plot.df <- plot.df[order(plot.df$log10Qval,decreasing = F),]
dim(plot.df)
if(nrow(plot.df)>20){
  plot.df <- tail(plot.df,20)
}

plot.df$Description <- factor(plot.df$Description,levels = plot.df$Description)
p <- ggplot(plot.df, aes(x = Description, y = log10Qval)) + 
  geom_bar(stat = 'identity', color = "darkgrey", fill = "darkgrey",width = 0.8) + 
  theme_classic() + coord_flip() +
  labs(x='',y = expression(paste("-log"[10], "(", italic("Q"), "-value)"))) +
  theme(axis.title = element_text(size = 8, colour = 'black'), 
        axis.text.y = element_text(size = 8, colour = 'black'), 
        axis.text.x = element_text(size = 8, colour = 'black'))
p
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/KEGG_',cutoff,'.pdf',sep=''),p,height = 2, width = 3)


### Tissue
Tissue <- read.table('enrichment.TISSUES.tsv',header = F,sep='\t')
plot.df <- Tissue
colnames(plot.df) <- c('BTO','Description','Counts','AllCounts','Strength','FDR','ENSG','Symbol')
print(plot.df$Description)
plot.df$log10Qval <- -log(plot.df$FDR,10)
plot.df <- plot.df[order(plot.df$log10Qval,decreasing = F),]
dim(plot.df)
if(nrow(plot.df)>10){
  plot.df <- tail(plot.df,10)
}

plot.df$Description <- factor(plot.df$Description,levels = plot.df$Description)
p <- ggplot(plot.df, aes(x = Description, y = log10Qval)) + 
  geom_bar(stat = 'identity', color = "darkgrey", fill = "darkgrey",width = 0.8) + 
  theme_classic() + coord_flip() +
  labs(x='',y = expression(paste("-log"[10], "(", italic("Q"), "-value)"))) +
  theme(axis.title = element_text(size = 8, colour = 'black'), 
        axis.text.y = element_text(size = 8, colour = 'black'), 
        axis.text.x = element_text(size = 8, colour = 'black'))
p
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/TISSUE_',cutoff,'.pdf',sep=''),p,height = 2, width = 3)


