rm(list=ls())
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggsignif)
library(xlsx)
library(reshape2)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)


plan='candiSet4'
cutoff='0'

# candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))

candiDB <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),
                      header  = T,sep='\t')
candiGenes <- unique(candiDB$gene)
length(candiGenes)


####################
#### TALL IC50
####################
GDSC <- read.csv('/local/yanzijun/CRU/TALL_FM/data/GDSC/GDSC2_fitted_dose_response_25Feb20.csv')
#GDSC <- read.csv('CRU/TALL_FM/data/GDSC/GDSC1_fitted_dose_response_25Feb20.csv')
print(length(unique(GDSC$COSMIC_ID)));print(length(unique(GDSC$DRUG_NAME)))
tmp <- GDSC[GDSC$TCGA_DESC=='ALL',c('COSMIC_ID','DRUG_NAME','LN_IC50')]
tmp$COSMIC_ID <- paste('DATA.',tmp$COSMIC_ID,sep='')
GDSC.ALL <- dcast(data=tmp,COSMIC_ID ~DRUG_NAME,fun.aggregate = sum)
rownames(GDSC.ALL) <- GDSC.ALL$COSMIC_ID;GDSC.ALL$COSMIC_ID=NULL
dim(GDSC.ALL)

####################
#### TALL celline exp
####################
meta <- read.xlsx('/local/yanzijun/CRU/TALL_FM/data/GDSC/Cell_Lines_Details.xlsx',sheetIndex = 1)
meta.TALL <- meta[meta$GDSC.Tissue.descriptor.2 %in% c('T_cell_leukemia','lymphoblastic_T_cell_leukaemia'),]
ID.TALL <- paste('DATA.',meta.TALL$COSMIC.identifier,sep='')
print(length(ID.TALL))

EXP <- read.csv('/local/yanzijun/CRU/TALL_FM/data/GDSC/Cell_line_RMA_proc_basalExp.txt',header = T,sep='\t') 
# EXP <- read.csv('/local/yanzijun/CRU/TALL_FM/data/DepMap/22Q1/CCLE_expression.csv')
# META <- read.csv('/local/yanzijun/CRU/TALL_FM/data/DepMap/21Q3/sample_info.csv')
dim(EXP)
EXP.bc <- EXP

EXP.TALL <- EXP[EXP$GENE_SYMBOLS %in% candiGenes,]
rownames(EXP.TALL) <- EXP.TALL$GENE_SYMBOLS
EXP.TALL$GENE_SYMBOLS=EXP.TALL$GENE_title=NULL
dim(EXP.TALL)
EXP.TALL <- as.data.frame(t(EXP.TALL[,colnames(EXP.TALL) %in% ID.TALL]))
EXP.TALL[1:3,1:4]
dim(EXP.TALL)


########################
### 1. 候选基因在TALL celline中的表达
########################
EXP.bg <- EXP[EXP$GENE_SYMBOLS %in% candiGenes,]
rownames(EXP.bg) <- EXP.bg$GENE_SYMBOLS
EXP.bg$GENE_SYMBOLS=EXP.bg$GENE_title=NULL
EXP.bg <- as.data.frame(t(EXP.bg[,!colnames(EXP.bg) %in% ID.TALL]))
EXP.bg[1:3,1:4]
dim(EXP.bg)
print(length(intersect(rownames(EXP.bg),rownames(EXP.TALL)))) #确认celline无交集


TALL.melt <- melt(EXP.TALL)
TALL.melt$Group='TALL'
bg.melt <- melt(EXP.bg)
bg.melt$Group='OtherCelline'
plot.df <- rbind(TALL.melt,bg.melt)
colnames(plot.df)[1:2] <- c('Gene','Exp')
plot.df$Gene=factor(plot.df$Gene,levels = candiGenes)
head(plot.df)

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/cellineEXP_',cutoff,'.pdf',sep=''),
    width = 10,height = 3)
p <- ggboxplot(plot.df, x = 'Group', y = "Exp",
               color = 'Group', 
               #palette = "jco",
               add = NULL,
               facet.by = "Gene", short.panel.labs = FALSE,nrow = 1)+
  stat_compare_means(label = "p.format")+
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust = 1 ,colour = 'black'))+
  scale_color_manual(values = getPalette[1:length(unique(plot.df$Group))])
print(p)
dev.off()




########################
### 2. 查看候选药物在T-ALL细胞系中的IC50分布
########################
Drug.df <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),header = T,sep='\t')
name2db <- Drug.df$DRUGBANK.ID
names(name2db) <- Drug.df$NAME
candiDrugs <- Drug.df$NAME

gene2drug <- Drug.df$gene
names(gene2drug) <- Drug.df$NAME


## 查看候选药物在GDSC中的交集，以及靶向的通路
Drug.ALL <- intersect(tmp$DRUG_NAME,candiDrugs)
print(Drug.ALL) #在ALL celline中的交集药物

drug2target <- read.csv('/local/yanzijun/CRU/TALL_FM/data/GDSC/screened_compunds_rel_8.2.csv',header = T)
df <- GDSC[GDSC$DRUG_NAME %in% Drug.ALL,
     c('DRUG_NAME','PUTATIVE_TARGET', 'PATHWAY_NAME')]
df <- df[!duplicated(df$DRUG_NAME),]
df$DrugBankID <- name2db[df$DRUG_NAME]
write.table(df,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/TargetPathway_',cutoff,'.txt',sep=''),
            col.names = T,row.names = F,sep='\t',quote = F)



candiGDSC <- tmp[tmp$DRUG_NAME %in%Drug.ALL ,]
candiGDSC$COSMIC_ID <- NULL
candiGDSC$Group='T-ALL'
head(candiGDSC)

## 背景值：在其他celline中，候选药物的IC50
bgGDSC <- GDSC[GDSC$DRUG_NAME %in% Drug.ALL & (!GDSC$COSMIC_ID %in% rownames(GDSC.ALL)),
               c('DRUG_NAME','LN_IC50')]
bgGDSC$Group='OtherCelline'
head(bgGDSC)

pvalue <- wilcox.test(exp(candiGDSC$LN_IC50),exp(bgGDSC$LN_IC50),alternative = 'less')
pvalue <- pvalue$p.value
print(pvalue)


#### plot
plot.df <- rbind(candiGDSC,bgGDSC)
plot.df$IC50 <- exp(plot.df$LN_IC50)
plot.df$DrugID <- name2db[plot.df$DRUG_NAME]
plot.df$Gene <- gene2drug[plot.df$DRUG_NAME]
head(plot.df)
plot.df <- plot.df[plot.df$Gene != 'WEE1',] ## 删除

plot.df$Tag <- paste(plot.df$DrugID,'\n(',plot.df$DRUG_NAME,')',sep='')
head(plot.df)

plot.df$Group=factor(plot.df$Group,levels = c("T-ALL","OtherCelline"))

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/IC50_',cutoff,'.pdf',sep=''),
    width = 3.5,height = 4)
p.ALL <- ggplot() +  geom_bar(data = plot.df, 
                              aes(x = Tag, y = IC50,fill=Group),
                              stat = "identity",
                              width = 0.7, 
                              position = position_dodge(width = 0.9)) + 
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'),
        axis.text = element_text(size=10,colour = 'black'),
        legend.position = 'bottom') + 
  #coord_flip()+
  scale_fill_manual(values = c(getPalette[1],'grey'))+
  #facet_wrap(~Gene,scales = "free")+
  labs(y="GDSC IC50",x='')
print(p.ALL)
dev.off()



### 截断
library(ggplot2)
library(ggpubr)
library("ggsci")
#画下面
p0 <- 
  ggplot(plot.df, aes(x = Tag, y = IC50, fill = Group))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  #stat_compare_means(method = "wilcox.test",method.args = list(alternative = "two.sided"))
  theme_classic2()+
  scale_fill_manual(values=c(getPalette[1],'grey'))+
    labs(x=NULL,y=NULL,fill=NULL)+    #可自定义标签名字
    coord_cartesian(ylim = c(0,250))+
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'),
          axis.text = element_text(size = 8,colour = 'black'),
          axis.title = element_text(size = 8,colour = 'black'),
          axis.line = element_line(size=0.5, colour = "black"),
          legend.position="none")+
    labs(y="GDSC IC50",x='')

#画上面
p1 <- ggplot(plot.df, aes(x = Tag, y = IC50, fill = Group))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))+
  labs(x=NULL,y=NULL,fill=NULL) +   #不要标签
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +     #去掉X轴和X轴的文字
  coord_cartesian(ylim = c(2000,5000)) +  #设置上面一半的值域
  scale_y_continuous(breaks = c(2000,3500,5000))+#以5为单位划分Y轴
  theme_classic2()+
  theme(axis.text.y = element_text(size = 8,colour = 'black'),
        text = element_text(size = 8,colour = 'black'),
        axis.line.y = element_line(size=0.5, colour = "black"),
        axis.line.x = element_line(size=0, colour = "white"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")+ 
  scale_fill_manual(values=c(getPalette[1],'grey'))+
  labs(x='',y='')

#拼起来
p <- ggarrange(p1,p0,heights=c(1/5, 4/5),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
p
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/IC50_',cutoff,'_cut.pdf',sep=''),p,
       width = 8,height = 7,units = 'cm')



########################
### 3. 候选基因在TALL celline中的表达与所有药物的IC50的相关性
########################
ol.Celline <- intersect(rownames(GDSC.ALL),rownames(EXP.TALL))
length(ol.Celline)

GDSC.df <- GDSC.ALL[match(ol.Celline,rownames(GDSC.ALL)),]
dim(GDSC.df)

EXP.df <- EXP.TALL[match(ol.Celline,rownames(EXP.TALL)),]
EXP.df <- EXP.df[match(rownames(EXP.df),rownames(GDSC.df)),]
print(all(rownames(GDSC.df)==rownames(EXP.df)))

cor_data <- cor(exp(GDSC.df),EXP.df,method = 'spearman')
dim(cor_data)
cor_data[1:3,1:5]
saveRDS(cor_data,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/GDSC/cor_',cutoff,'.RDS',sep=''))




saveRDS(cor_data,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/GDSC/cor_',cutoff,'.RDS',sep=''))

#### plot 候选药物的GDSC相关性，没什么意义，都在各自的靶点里高相关性
cor.ALL <- cor_data[Drug.ALL,]
p.GDSC <- pheatmap::pheatmap(cor.ALL,cluster_cols = T,cluster_rows = T, 
                             show_colnames = T,
                             border_color = "white",
                             colorRampPalette(c("navy", "white", "firebrick3"))(50),
)

ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/corGDSC_',cutoff,'.pdf',sep=''),p.GDSC,height = 25,width = 6)


#########################
## GDSC
#########################
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(ggcorrplot)
library(pheatmap)

## 过滤药物
func_filter <- function(df,cut_p=0.5,pct=7){
  index <- c()
  for(i in 1:nrow(df)){
    x=df[i,]
    n1 <- length(x[x > cut_p])
    n2 <- length(x[x < -cut_p])
    if(n1>length(x)/pct | n2>length(x)/pct){
      index <- c(index,i)
    }
  }
  print(length(index))
  
  new.df <- df[index,]
  return(new.df)
}

drugs <- list()
for(i in 1:ncol(cor_data)){
  d <- cor_data[,i]
  # d <- sort(d,decreasing = T)
  # d <- c(head(d,10),tail(d,10))
  d <- d[abs(d)>0.4]
  drugs[[i]] <- names(d)
}
freq <- table(unlist(drugs))
candiDrugs <- names(freq[as.numeric(freq)>=3])
print(length(candiDrugs))

new.df <- cor_data[rownames(cor_data) %in% candiDrugs,]
dim(new.df)

p.GDSC <- pheatmap::pheatmap(cor.ALL,cluster_cols = T,cluster_rows = T, 
                             show_colnames = T,
                             border_color = "white",
                             colorRampPalette(c("navy", "white", "firebrick3"))(50),
)


