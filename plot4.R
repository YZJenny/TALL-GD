rm(list=ls())

library(foreach)
library(reshape2)
library(tibble)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
coul  <- brewer.pal(3, "Set1")

plan='candiSet4'
cutoff='0'

candiDB <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),
            header  = T,sep='\t')
candiDB$FDA <- 'Unapproved'
candiDB$FDA[grep('approved',candiDB$DRUG.GROUP)] <- 'Approved'

#link <- candiDB[,c('Gene','DRUGBANK.ID')]
link <- candiDB[,c('DRUGBANK.ID','gene')]
nodes.Drug <- data.frame(name=candiDB$DRUGBANK.ID,type=candiDB$FDA,group=rep('Drug',nrow(candiDB)))
nodes.Gene <- data.frame(name=candiDB$gene,type=rep('AGene',nrow(candiDB)),group=rep('Gene',nrow(candiDB)))
nodes <- rbind(nodes.Drug,nodes.Gene)
nodes <- nodes[!duplicated(nodes$name),]

g <- graph_from_data_frame(d = link,vertices = nodes,directed = TRUE)
my_color <- coul[as.numeric(as.factor(V(g)$type))]
my_shape <- c("square",'sphere')[as.numeric(as.factor(V(g)$group))]

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/network_',cutoff,'.pdf',sep='')
    ,width = 10,height = 10)
# 绘值带节点属性的网络图
plot(g, 
     vertex.color=my_color,
     vertex.frame.color='white',
     vertex.shape=my_shape,
     vertex.size=3,
     vertex.label.color='black',
     vertex.label.dist=0.5,
     vertex.label.degree=30,
     edge.color='grey',
     edge.width=1,           
     edge.arrow.size=0.5,           
     edge.arrow.width=1,                
     edge.lty=c("solid"))
# 添加图例
legend("bottomleft", legend=levels(as.factor(V(g)$type)), 
       col = coul , bty = "n", pch=20 , pt.cex = 3, 
       cex = 1, text.col=coul , horiz = FALSE, 
       inset = c(0, 0))
dev.off()
