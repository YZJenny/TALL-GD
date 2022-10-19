############  
### Clusterprofiler
############
rm(list=ls())
#.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",'/local/yanzijun/tools/R-4.0.0/library'))
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(AnnotationHub)
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(DOSE)
library(dplyr)
library(tidyverse)
library(reshape2)

plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

geneLst <- candiGenes
symbol2id=mapIds(org.Hs.eg.db,geneLst,"ENTREZID",'SYMBOL')
id=symbol2id[which(symbol2id!='')] #提取出非NA的ENTREZID
print(length(id))
#GO富集分析#
ego <- enrichGO(OrgDb="org.Hs.eg.db", gene = id, ont = "BP", 
                pvalueCutoff = 0.05, readable= TRUE) #GO富集分析
print(dim(ego))
saveRDS(ego,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/GO_',cutoff,'_all.RDS',sep=''))

ego.simp <- clusterProfiler::simplify(ego,cutoff=0.6,by="p.adjust",select_fun=min,measure="Wang")
print(dim(ego.simp))
ego_res <- as.data.frame(ego.simp)
print(ego_res$Description)
saveRDS(ego.simp,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/GO_',cutoff,'.RDS',sep=''))
write.table(ego_res,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/GO_',cutoff,'.txt',sep=''),sep='\t',quote=F,row.names = F)


#KEGG分析#
ekk <- enrichKEGG(gene= id,organism  = 'hsa',pvalueCutoff=0.05)	 #KEGG富集分析
print(dim(ekk))
if(nrow(ekk)>0){
  saveRDS(ekk,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/KEGG_',cutoff,'.RDS',sep=''))
  ekk_res <- as.data.frame(ekk)
  print(ekk_res$Description)
  write.table(ekk_res,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/KEGG_',cutoff,'.txt',sep=''),sep='\t',quote=F,row.names = F)
}


########
## 画图
########
library(xlsx)

ego_res <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Enrich/GO_',cutoff,'.RDS',sep=''))
ego_res <- as.data.frame(ego_res)
plot.df <- ego_res
#plot.df <- ekk_res
print(plot.df$Description)

plot.df$log10Qval <- -log(plot.df$qvalue,10)
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
  theme(axis.title = element_text(size = 10, colour = 'black'), 
        axis.text.y = element_text(size = 10, colour = 'black'), 
        axis.text.x = element_text(size = 10, colour = 'black'))
p
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/GO_',cutoff,'.pdf',sep=''),p,height = 6, width = 6)
#ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/KEGG_',cutoff,'.pdf',sep=''),p,height = 3, width = 5)

