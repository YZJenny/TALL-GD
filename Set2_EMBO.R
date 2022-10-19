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
## EMBO 筛选1/3象限的基因和选差异的ATAC peak
############
## DAP和DEG，筛选1/3象限的基因
all.DAEG <- read.csv('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/allDAEG_EV9.csv')
all.DAEG <- all.DAEG[1:4163,1:5]
head(all.DAEG)

table(all.DAEG$quadrant)
up.Genes <- all.DAEG$gene[all.DAEG$quadrant==1]
dn.Genes <- all.DAEG$gene[all.DAEG$quadrant==3]

print(length(up.Genes))
print(length(dn.Genes))

GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'CRU/TALL_FM/res/GeneSet/EMBO/GeneLst_paper.RDS')


#### plot
fc_cutoff=1
data <- all.DAEG
colnames(data)[2:3] <- c('log2FC_RNA','log2FC_ATAC')
data$quadrant <- factor(data$quadrant)
figure <- ggplot(data,aes(x=log2FC_RNA,y=log2FC_ATAC))+
  geom_point(aes(color=quadrant))+
  scale_color_manual(values = c("#FF0000","#999999","#999999","#999999"))+
  theme_bw()+
  labs(title='Busra et al',x="log2 (Fold Change) RNA",y = "log2 (Fold Change) ATAC")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12,colour = 'black'),
        axis.line = element_line(size=0.5),
        panel.grid.major=element_line(colour=NA),panel.grid=element_blank())
p <- figure+
  geom_vline(xintercept=log(fc_cutoff,2),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = log(fc_cutoff,2),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
p
ggsave('/local/yanzijun/CRU/TALL_FM/res/FIG/Quadrant_paper.pdf',p)


###############
### 作者提供的文件，选差异的ATAC peak
###############
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

all.peak <- read.csv('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/DAP_EV11.csv')

up.peak <- all.peak[all.peak$log2FoldChange>0,]
dn.peak <- all.peak[all.peak$log2FoldChange<0,]
print(dim(up.peak))
print(dim(dn.peak))

up.Genes <- c()
for(i in 1:nrow(up.peak)){
  symbols <- unlist(strsplit(up.peak$annotation[i],','))
  up.Genes <- c(up.Genes,symbols)
}

dn.Genes <- c()
for(i in 1:nrow(dn.peak)){
  symbols <- unlist(strsplit(dn.peak$annotation[i],','))
  dn.Genes <- c(dn.Genes,symbols)
}

print(length(up.Genes))
print(length(dn.Genes))

GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/GeneLst.RDS')



#### plot
data <- all.peak
data <- data[,-c(ncol(data)-1,ncol(data))]
qvalue <- 0.05
fc_cutoff <- 1
data <- mutate(data,
               sig=ifelse((data$padj < qvalue & data$log2FoldChange > log(fc_cutoff,2))|
                            (data$padj < qvalue & data$log2FoldChange < -log(fc_cutoff,2)) ,
                          ifelse(data$log2FoldChange > log(fc_cutoff,2),'UP','DOWN'),'no'))
data <- na.omit(data)
data$log10Qvalue=-log10(data$padj)

figure <- ggplot(data,aes(x=log2FoldChange,y=log10Qvalue))+geom_point(aes(color=sig))+
  scale_color_manual(values = c(mycol[2],mycol[1]))+
  theme_bw()+
  labs(title='Busra et al',x="log2 (Fold Change)",y = "-log10 (adjP value)")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12,colour = 'black'),
        axis.line = element_line(size=0.5),
        panel.grid.major=element_line(colour=NA),panel.grid=element_blank())
p <- figure+
  geom_vline(xintercept=log(fc_cutoff,2),lty=3,col="black",lwd=0.5) +#添加横线|FoldChange|>2
  geom_hline(yintercept = -log10(qvalue),lty=3,col="black",lwd=0.5)#添加竖线padj<0.05
p
ggsave('/local/yanzijun/CRU/TALL_FM/res/FIG/DAP_EV11.pdf',p,width = 4,height = 3.5)



##############################
### 自己用Limma算(EV7)
##############################
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0",
            "/usr/local/lib64/R/library"))
library(limma)

peak.df <- read.csv('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/peak_EV7.csv',header  = T)
colnames(peak.df)[7:9]

data <- peak.df[,8:ncol(peak.df)]
rownames(data) <- peak.df$id
label <- rep('tumor',ncol(data))
label[grep('^T',colnames(data))]='normal'
Group <- data.frame(ID=colnames(data),ID_REF=label)

write.table(data,'/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/Peak_Value.txt',quote = F,sep='\t',row.names = T,col.names = T)
write.table(Group,'/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/Peak_design.txt',quote = F,sep='\t',row.names = F,col.names = F)


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
output <- tibble::rownames_to_column(output,'id')  
dim(output)
output <- merge(peak.df[,c('id','annotation_geneSymbol')],output,by='id')
write.table(output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/Limma_EV7.txt',sep='\t',row.names = F,quote = F)


### 
output <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/Limma_EV7.txt',sep='\t',header = T)
peak.df <- read.csv('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/peak_EV7.csv',header  = T)
print(all(output$id==peak.df$id))
output$annotation_geneSymbol<- as.character(output$annotation_geneSymbol)
output$peak_type <- peak.df$peak.type..TSS.or.distal.
output$peak_category <- peak.df$peak.category

fc_cutoff=1.3
up.peak <- output[output$logFC>log(fc_cutoff,2) & output$adj.P.Val < 0.05,]
dn.peak <- output[output$logFC< -log(fc_cutoff,2) & output$adj.P.Val < 0.05,]
print(dim(up.peak))
print(dim(dn.peak))

up.Genes <- c()
for(i in 1:nrow(up.peak)){
  print(i)
  symbols <- unlist(strsplit(up.peak$annotation_geneSymbol[i],','))
  up.Genes <- unique(c(up.Genes,symbols))
}

dn.Genes <- c()
for(i in 1:nrow(dn.peak)){
  print(i)
  symbols <- unlist(strsplit(dn.peak$annotation_geneSymbol[i],','))
  dn.Genes <- unique(c(dn.Genes,symbols))
}

print(length(up.Genes))
print(length(dn.Genes))

GeneLst <- list(up.Genes=up.Genes,dn.Genes=dn.Genes)
saveRDS(GeneLst,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/GeneLst_EV7.RDS')


### 对应基因的limma res
gene_output <- c()
for(i in 1:nrow(output)){
  print(i)
  symbols <- unlist(strsplit(output$annotation_geneSymbol[i],','))
  for(j in 1:length(symbols)){
    g=symbols[j]
    df <- data.frame(gene=g,logFC=output$logFC[i])
    gene_output <- rbind(gene_output,df)
  }
}
saveRDS(gene_output,'/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/Limma_EV7_gene.RDS')


##############
## Heatmap
##############
rm(list=ls())
library(ggplot2)
library(pheatmap)
limma_res <- read.table('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/Limma_EV7.txt',sep='\t',header = T)
peak.df <- read.csv('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/peak_EV7.csv',header  = T)
limma_res$peak.type <- peak.df$peak.type..TSS.or.distal.


EXP <- read.table('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/Peak_Value.txt',header=T,sep='\t')
Group <- read.table('/local/yanzijun/CRU/TALL/data/EMBO_2020/suppment/Peak_design.txt',header=F,sep='\t')
Label=Group[,2]

FC=1.3
up.sigG <- limma_res[which(limma_res$adj.P.Val < 0.05 & limma_res$logFC > log(FC,2) ),]
dim(up.sigG)
top.up.sigG <- up.sigG$id[order(up.sigG$logFC,decreasing = T)][1:100]

dn.sigG <- limma_res[which(limma_res$adj.P.Val < 0.05 & limma_res$logFC < -log(FC,2) ),]
dim(dn.sigG)
top.dn.sigG <- dn.sigG$id[order(dn.sigG$logFC,decreasing =F )][1:100]
top.sigG <- c(top.up.sigG,top.dn.sigG)
sigG <- c(up.sigG$id,dn.sigG$id)

top.sigG <- sigG
#sigG <- candiGenes
print(length(top.sigG))

## heatmap
df_heatmap <- EXP[which(rownames(EXP)%in% top.sigG),]
print(dim(df_heatmap))

annotation_col = data.frame(SampleType = factor(Label))
rownames(annotation_col) = colnames(df_heatmap)


plot.exp <- df_heatmap[match(top.sigG,rownames(df_heatmap)),] #行匹配
print(all(rownames(plot.exp)==top.sigG))


ann_colors = list(
  SampleType=c(normal=mycol[2],tumor=mycol[1]))

p.heatmap <- pheatmap::pheatmap(plot.exp,scale = 'row',cluster_cols = T,cluster_rows = F, 
                                annotation_col = annotation_col,show_colnames = F,show_rownames = F,
                                treeheight_row=0,treeheight_col=12,
                                annotation_colors = ann_colors,border_color = "white",
                                colorRampPalette(c("navy","white",  "firebrick3"))(length(seq(-3,3,by = 0.1))),
                                gaps_col = cumsum(
                                  c(length(which(Label=='normal')),length(which(Label=='tumor'))) ),
                                breaks = seq(-3,3,by = 0.1),legend_breaks = seq(-3,3,1))

ggsave('/local/yanzijun/CRU/TALL_FM/res/candiSet4/FIG/DAP_heatmap.pdf',p.heatmap,width = 4,height = 6)


##################
### barplot
##################
Freq.num <- data.frame(type=c('open','close'),
                       Freq=c(nrow(up.sigG),-nrow(dn.sigG)))
Freq.num
p.Freq <- ggplot(data=Freq.num, mapping=aes(x=type,y=Freq))+
  geom_bar(stat="identity")+
  theme_classic()
ggsave('/local/yanzijun/CRU/TALL_FM/res/candiSet4/FIG/DAP_Freq.pdf',p.Freq,width = 2,height = 4)



####
sig.df <- limma_res[limma_res$id %in% sigG,]
sig.df[1:2,]
sig.df$type <- 'open'
sig.df$type[sig.df$id %in% dn.sigG$id]='close'
sig.df <- sig.df[,c('peak.type','type')]
table(sig.df$peak.type,sig.df$type)

mytable <- with(sig.df,table(peak.type,type))
df <- as.data.frame(mytable)
df$Freq[1:2] <- -df$Freq[1:2]
df

p.Num <- ggplot(data=df, mapping=aes(x=peak.type,y=Freq,fill=type))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  theme_classic()+
  scale_fill_manual(values= c("#999999",mycol[1]))
p.Num
ggsave('/local/yanzijun/CRU/TALL_FM/res/candiSet4/FIG/DAP_Num_FC1.3.pdf',p.Num,width = 5,height = 4)
