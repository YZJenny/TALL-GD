##### 候选基因在临床信息的显著性
rm(list=ls())
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(pheatmap)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)


## gene exp file
EXP <- as.data.frame(read_excel("/local/yanzijun/CRU/TALL/data/NG_2017/RNAseq_FPKM.xlsx"))
rownames(EXP) <- EXP$Gene;EXP$Gene <- NULL
EXP <- log(EXP+0.1,2)

## clinical file
raw_info <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/TARGET_NG2017.txt',header = T,sep='\t')

ol.samples <- intersect(colnames(EXP),raw_info$RNAseq_id_D)
raw_info <- raw_info[raw_info$RNAseq_id_D %in% ol.samples,c('RNAseq_id_D','Age_Dx_yrs','Gender',
                                                              'WBC.at.Diagnosis','CNS.Status.at.Diagnosis',
                                                              'Bone.Marrow.Site.of.Relapse','CNS.Site.of.Relapse',
                                                            'ETP.status','BM_blasts_diagnosis','group')]
rownames(raw_info) <- raw_info$RNAseq_id_D

EXP <- select(EXP,raw_info$RNAseq_id_D)
print(all(colnames(EXP)==raw_info$RNAseq_id_D))


new_info <- data.frame(sampleID=raw_info$RNAseq_id_D)

new_info$Age <- as.numeric(raw_info$Age_Dx_yrs)
print(summary(new_info$Age))

new_info$Gender=raw_info$Gender
# new_info$Gender <- 1
# new_info[which(raw_info$Gender=='Female'),ncol(new_info)] <- 0

new_info$WBC <- 'WBC1'
new_info[as.numeric(raw_info$WBC.at.Diagnosis)>50 & as.numeric(raw_info$WBC.at.Diagnosis) <200,ncol(new_info)] <- 'WBC2'
new_info[which(as.numeric(raw_info$WBC.at.Diagnosis)>200),ncol(new_info)] <- 'WBC3'
new_info$WBC <- factor(new_info$WBC,levels = c('WBC1','WBC2','WBC3'))

# new_info$WBC <- 'High'
# new_info[as.numeric(raw_info$WBC.at.Diagnosis)< median(as.numeric(raw_info$WBC.at.Diagnosis)),ncol(new_info)] <- 'Low'
# new_info$WBC <- factor(new_info$WBC,levels = c('Low','High'))


new_info$CNS.status <- 'CNS1'
new_info[which(raw_info$CNS.Status.at.Diagnosis %in% c('CNS 2a','CNS 2b','CNS 2c')),ncol(new_info)] <- 'CNS2'
new_info[which(raw_info$CNS.Status.at.Diagnosis %in% c('CNS 3a','CNS 3b','CNS 3c')),ncol(new_info)] <- 'CNS3'
new_info$CNS.status <- factor(new_info$CNS.status,levels = c('CNS1','CNS2','CNS3'))


new_info$BM.Relapse <- 'Yes'
new_info[which(raw_info$Bone.Marrow.Site.of.Relapse=='No'),ncol(new_info)] <- 'No'
new_info$BM.Relapse <- factor(new_info$BM.Relapse,levels = c('Yes','No'))


new_info$CNS.Relapse <- 'Yes'
new_info[which(raw_info$CNS.Site.of.Relapse=='No'),ncol(new_info)] <- 'No'
new_info$CNS.Relapse <- factor(new_info$CNS.Relapse,levels = c('Yes','No'))

new_info$ETP.status <- NA
new_info[which(raw_info$ETP.status=='ETP'),ncol(new_info)] <- 'ETP'
new_info[which(raw_info$ETP.status=='nearETP'),ncol(new_info)] <- 'nearETP'
new_info[which(raw_info$ETP.status=='notETP'),ncol(new_info)] <- 'notETP'
table(new_info$ETP.status)
new_info$ETP.status <- factor(new_info$ETP.status,levels = c('ETP','nearETP','notETP'))


new_info$BM_blasts_diagnosis <- 'High'
new_info[which(raw_info$BM_blasts_diagnosis < median(raw_info$BM_blasts_diagnosis)),ncol(new_info)] <- 'Low'
table(new_info$BM_blasts_diagnosis)
new_info$BM_blasts_diagnosis <- factor(new_info$BM_blasts_diagnosis,levels = c('High','Low'))

new_info$Age_b <- 'children'
new_info$Age_b[new_info$Age>12]='Adult'

## add subtype info
meta <- as.data.frame(read_excel("/local/yanzijun/CRU/TALL/data/NG_2017/cohortInfo.xlsx"))
table(meta$group)
subtype <- meta$group
names(subtype) <- meta$RNAseq_id_D
new_info$subtype= subtype[new_info$sampleID]
table(new_info$subtype)
dim(new_info)
new_info[1:2,]

### 
plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)


candi.EXP <- as.data.frame(t(EXP))
candi.EXP <- candi.EXP[,colnames(candi.EXP)%in%candiGenes]

print(all(rownames(candi.EXP)==new_info$sampleID))


################
#### boxplot 
################
Variables <- c('ETP.status','subtype','BM_blasts_diagnosis','WBC','CNS.status','CNS.Relapse','BM.Relapse','Age_b')
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/Clinical_',cutoff,'.pdf',sep=''),
    width = 12,height = 7)
for(i in 1:length(Variables)){
  Var=Variables[i]
  print(Var)
  df <- candi.EXP
  df[,Var] <- new_info[,Var]
  df <- na.omit(df)
  
  plot.df <- reshape2::melt(df,value.name = Var)
  head(plot.df)
  colnames(plot.df)[2:3] <- c('Gene','Exp')
  
  p <- ggboxplot(plot.df, x = Var, y = "Exp",
                 color = Var, 
                 #palette = "jco",
                 add = NULL,
                 facet.by = "Gene", short.panel.labs = FALSE,nrow = 3)+
    stat_compare_means(label = "p.signif",
                       label.x = 4,label.y=9)+
    scale_color_manual(values = mycol[1:length(unique(plot.df[,Var]))])
  print(p)
  
}
dev.off()
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/TARGETsubtype_',cutoff,'_boxplot.pdf',sep=''),p,
       width = 11,height = 7)
    

### 只选取 ETP显著的展示
Var='ETP.status'
print(Var)

sig.Genes <- c('ACACA','MPZL1','STT3B','YWHAG','GNAQ','BUD13','CUEDC2','DNAJC9','HMGCS1','POLD2','PRMT5','PRPF19',
               'PSMB2','RBBP8','RCN1','SCD','TUBB','TYMS','ZNF273')

## significant genes
df <- candi.EXP[,colnames(candi.EXP)%in%sig.Genes]
df[,Var] <- new_info[,Var]
df <- na.omit(df)

plot.df <- melt(df,value.name = Var)
head(plot.df)
colnames(plot.df)[2:3] <- c('Gene','Exp')
plot.df$Gene <- factor(plot.df$Gene,levels = sig.Genes)

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/Clinical_',cutoff,'_ETP_sig.pdf',sep=''),
    width = 7,height = 9)
p <- ggboxplot(plot.df, x = Var, y = "Exp",
               color = Var, 
               #palette = "jco",
               #add = "jitter",
               add=NULL,
               facet.by = "Gene", short.panel.labs = FALSE,nrow = 4)+
  stat_compare_means(label = "p.format")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'))+
  scale_color_manual(values = mycol[1:length(unique(plot.df[,Var]))])
print(p)
dev.off()




## NOT significant genes
df <- candi.EXP[,!colnames(candi.EXP)%in%sig.Genes]
df[,Var] <- new_info[,Var]
df <- na.omit(df)

plot.df <- melt(df,value.name = Var)
head(plot.df)
colnames(plot.df)[2:3] <- c('Gene','Exp')

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/Clinical_',cutoff,'_ETP_NOTsig.pdf',sep=''),
    width = 7,height = 5.5)
p <- ggboxplot(plot.df, x = Var, y = "Exp",
               color = Var, 
               #palette = "jco",
               add = NULL,
               facet.by = "Gene", short.panel.labs = FALSE,nrow = 2)+
  stat_compare_means(label = "p.format")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1 ,colour = 'black'))+
  scale_color_manual(values = mycol[1:length(unique(plot.df[,Var]))])
print(p)
dev.off()


##################
## heatmap:ETP
##################
ETP <- new_info$sampleID[which(new_info$ETP.status=='ETP')]
nearETP <- new_info$sampleID[which(new_info$ETP.status=='nearETP')]
notETP <- new_info$sampleID[which(new_info$ETP.status=='notETP')]

sample_order <- c(ETP,nearETP,notETP)
plot.exp <- dplyr::select(as.data.frame(t(candi.EXP)),sample_order) #列匹配
print(all(colnames(plot.exp)==sample_order))


library(RColorBrewer)
library(ggplot2)
col2 =brewer.pal(9,"Set1")[1:9]

group_df = data.frame(
  Subtype=rep(c('ETP','nearETP','notETP' ), 
              c(length(ETP),length(nearETP),length(notETP))))
rownames(group_df) <- colnames(plot.exp)
ann_colors = list(
  Subtype=c(ETP=col2[1],nearETP=col2[2],notETP=col2[3]))

p.heatmap <- pheatmap::pheatmap(plot.exp,scale='row',cluster_cols = F,cluster_rows = T, 
                               annotation_col = group_df,show_colnames = F,
                               treeheight_row=0,treeheight_col=0,
                               annotation_colors = ann_colors,border_color = "white",
                               colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(-4,4,by = 0.1))),
                               gaps_col = cumsum(
                                 c(length(ETP),length(nearETP),length(notETP))),
                               breaks = seq(-4,4,by = 0.1),legend_breaks = seq(-4,4,2))

ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/ETP_subtype_',cutoff,'.pdf',sep=''),p.heatmap,width = 6,height = 6)



##################
## heatmap:subtype
##################
TAL1 <- new_info$sampleID[new_info$subtype=='TAL1']
TAL2 <- new_info$sampleID[new_info$subtype=='TAL2']
HOXA <- new_info$sampleID[new_info$subtype=='HOXA']
LMO1_2 <- new_info$sampleID[new_info$subtype=='LMO1/2']
LMO2_LYL1 <- new_info$sampleID[new_info$subtype=='LMO2_LYL1']
NKX2_1 <- new_info$sampleID[new_info$subtype=='NKX2_1']
TLX1 <- new_info$sampleID[new_info$subtype=='TLX1']
TLX3 <- new_info$sampleID[new_info$subtype=='TLX3']
Unknown <- new_info$sampleID[new_info$subtype=='Unknown']

sample_order <- c(TAL1,TAL2,HOXA,LMO1_2,LMO2_LYL1,NKX2_1,TLX1,TLX3,Unknown )
plot.exp <- dplyr::select(as.data.frame(t(candi.EXP)),sample_order) #列匹配
print(all(colnames(plot.exp)==sample_order))


library(RColorBrewer)
library(ggplot2)
col2 =brewer.pal(9,"Set1")[1:9]

group_df = data.frame(
  Subtype=rep(c('TAL1','TAL2','HOXA','LMO1_2','LMO2_LYL1','NKX2_1','TLX1','TLX3','Unknown' ), 
              c(length(TAL1),length(TAL2),length(HOXA),
                length(LMO1_2),length(LMO2_LYL1),length(NKX2_1),length(TLX1),length(TLX3),length(Unknown))))
rownames(group_df) <- colnames(plot.exp)
ann_colors = list(
  Subtype=c(TAL1=col2[1],TAL2=col2[2],HOXA=col2[3],LMO1_2=col2[4],LMO2_LYL1=col2[5],
            NKX2_1=col2[6],TLX1=col2[7],TLX3=col2[8],Unknown=col2[9]))


p.TARGET <- pheatmap::pheatmap(plot.exp,scale = 'row',cluster_cols = F,cluster_rows = T, 
                               annotation_col = group_df,show_colnames = F,
                               treeheight_row=0,treeheight_col=0,
                               annotation_colors = ann_colors,border_color = "white",
                               colorRampPalette(c("navy", "white", "firebrick3"))(length(seq(-4,4,by = 0.1))),
                               gaps_col = cumsum(
                                 c(length(TAL1),length(TAL2),length(HOXA),
                                   length(LMO1_2),length(LMO2_LYL1),length(NKX2_1),length(TLX1),length(TLX3),length(Unknown))
                               ),
                               breaks = seq(-4,4,by = 0.1),legend_breaks = seq(-4,4,2))
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/TARGETsubtype_',cutoff,'_heatmap.pdf',sep=''),p.TARGET,width = 10,height = 4)



##################
## heatmap:WBC
##################
WBC1 <- new_info$sampleID[which(new_info$WBC=='WBC1')]
WBC2 <- new_info$sampleID[which(new_info$WBC=='WBC2')]
WBC3 <- new_info$sampleID[which(new_info$WBC=='WBC3')]

sample_order <- c(WBC1,WBC2,WBC3)
plot.exp <- dplyr::select(as.data.frame(t(candi.EXP)),sample_order) #列匹配
print(all(colnames(plot.exp)==sample_order))


library(RColorBrewer)
library(ggplot2)
col2 =brewer.pal(9,"Set1")[1:9]

group_df = data.frame(
  Subtype=rep(c('WBC1','WBC2','WBC3' ), 
              c(length(WBC1),length(WBC2),length(WBC3))))
rownames(group_df) <- colnames(plot.exp)
ann_colors = list(
  Subtype=c(WBC1=col2[1],WBC2=col2[2],WBC3=col2[3]))


p.heatmap <- pheatmap::pheatmap(plot.exp,scale = 'row',cluster_cols = F,cluster_rows = T, 
                                annotation_col = group_df,show_colnames = F,
                                treeheight_row=0,treeheight_col=0,
                                annotation_colors = ann_colors,border_color = "white",
                                colorRampPalette(c("navy","white",  "firebrick3"))(length(seq(-4,4,by = 0.1))),
                                gaps_col = cumsum(
                                  c(length(WBC1),length(WBC2),length(WBC3)) ),
                                breaks = seq(-4,4,by = 0.1),legend_breaks = seq(-4,4,2))

ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/WBC_subtype_',cutoff,'.pdf',sep=''),p.heatmap,width = 10,height = 4)


##########
## survival
##########
rm(list=ls())
.libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
            "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
library(ggplot2)
library(dplyr)
library("survival")
library("survminer")
library(survivalROC)
library(rmda)
library(DCA)
library(RColorBrewer)
mycol = colorRampPalette(brewer.pal(9, "Set1"))(9)

### 
plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

EXP <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/inputExp.txt',header = T,sep='\t')
EFS <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/EFS_subtype.txt',header = T,row.names = 1,sep='\t')
OS <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/OS_subtype.txt',header = T,row.names = 1,sep='\t')
Surv <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/GeneLst.RDS')
print(length(Surv))


group='sig'
if(group=='sig'){
  Genes <- intersect(candiGenes,Surv)
}else if(group=='NOTsig'){
  Genes <- setdiff(candiGenes,Surv)
}


Type=c('EFS','OS')
for(j in 1:length(Type)){
  type=Type[j]
  if(type=='EFS'){
    CLIN=EFS
  }else if(type=='OS'){
    CLIN=OS
  }
  CLIN$surtype <- NULL
  
  EXP_tmp <- as.data.frame(t(EXP))
  EXP_tmp <- EXP_tmp[,sapply(EXP_tmp, function(x)
    ifelse(sum(x != 0, na.rm = TRUE) > 30, TRUE, FALSE)) ] #删除列中包含>30个0的基因，否则会报错
  EXP_tmp <- tibble::rownames_to_column(EXP_tmp,'ID')
  CLIN_tmp <- tibble::rownames_to_column(CLIN,'ID')
  data <- merge(CLIN_tmp,EXP_tmp)
  
  res.cut <- surv_cutpoint(data,time = "time", event = "status",variables = Genes)
  res.cat <- surv_categorize(res.cut)
  write.table(summary(res.cut),
              paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/Cut_',type,'_',group,'.txt',sep=''),
              row.names = T,col.names = T,quote = F,sep='\t')
  
  write.table(res.cat,
              paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/Cat_',type,'_',group,'.txt',sep=''),
              row.names = F,col.names = T,quote = F,sep='\t')
  
  clin_info=res.cat[,1:2]
  group_info <- res.cat[,3:ncol(res.cat)]
  Genes=colnames(group_info)
  
  pvalue.plots <- list()
  for(k in 1:length(Genes)){
    surtype=group_info[,k]
    clin_info$surtype <- surtype
    fit <- survfit(Surv(time =time, event =status) ~ surtype,data=clin_info)
    #pvalue <- surv_pvalue(fit)$pval
    
    pvalue.plots[k] <- ggsurvplot(fit,data = clin_info, 
                                  #risk.table = TRUE,
                                  pval = T,palette = 'nejm',
                                  title = paste(type,' of ',Genes[k],sep=''),
                                  ggtheme = theme_survminer(size = 1))
    # pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/',Genes[k],'_',type,'_',cutoff,'.pdf',sep=''),
    #     width = 6,height = 6)
    # print(p, newpage = FALSE)
    # dev.off()
  }
  pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/',type,'_',group,'.pdf',sep=''),width = 3,height = 3)
  print(pvalue.plots)
  dev.off()
  
  p <- ggarrange(plotlist = pvalue.plots,nrow = 4,ncol=2)
  print(p)
  ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/',type,'_',group,'_2.pdf',sep=''),p,width = 6,height = 13)
}
