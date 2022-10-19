rm(list=ls())

library(xlsx)
library(ggplot2)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)
show_col(getPalette)

## candi 5 target Genes
plan='candiSet4'
cutoff='0'

candiDrug.df <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),header = T,sep='\t')
geneLst <- c('CDK6','TUBB','TUBA1A','TYMS','PSMB2')
candiDrug.df <- candiDrug.df[candiDrug.df$gene %in% geneLst,]
dim(candiDrug.df)

candiDrug <- unique(candiDrug.df$NAME)
print(length(candiDrug)) # 候选药物个数
approved=candiDrug.df[grep('approve',candiDrug.df$DRUG.GROUP),]
print(length(unique(approved$DRUGBANK.ID))) # 候选药物个数

unapproved=candiDrug.df[-grep('approve',candiDrug.df$DRUG.GROUP),]
print(length(unique(unapproved$DRUGBANK.ID))) # 候选药物个数


candiDrugID <- candiDrug.df$DRUGBANK.ID
names(candiDrugID) <- candiDrug.df$NAME

candiGene <- candiDrug.df$gene
names(candiGene) <- candiDrug.df$NAME


## NCI60  
NCI60.score <- read.csv('/local/yanzijun/public/NCI60/DTP_NCI60_ZSCORE.csv',header = T)
dim(NCI60.score)
NCI60.score[1:2,]
NCI60.score$NSC=NCI60.score$FDA.status=NCI60.score$Mechanism.of.action.c=
  NCI60.score$SMILES.d=NCI60.score$Total.experiments.e=NCI60.score$Total.after.quality.control.f=NULL
NCI60.score[1:2,]


# ## bg: other drugs
# bg.score <- NCI60.score[-which(NCI60.score$Drug.name %in% candiDrug | 
#                             NCI60.score$PubChem.SID %in% as.character(candiDrug.df$CID)),
#                           c('Drug.name','LE.CCRF.CEM','LE.MOLT.4')]
# 
# dim(bg.score)
# print(length(unique(bg.score$Drug.name))) #有重复实验，算均值
# bg.score$LE.CCRF.CEM <- as.numeric(bg.score$LE.CCRF.CEM)
# bg.score$LE.MOLT.4 <- as.numeric(bg.score$LE.MOLT.4)
# dim(bg.score)
# 
# bg.CCRF_CEM<- bg.score[,c('Drug.name','LE.CCRF.CEM')]
# bg.CCRF_CEM <- na.omit(bg.CCRF_CEM)
# bg.CCRF_CEM <-   aggregate(bg.CCRF_CEM$LE.CCRF.CEM, by=list(type=bg.CCRF_CEM$Drug.name),mean)
# colnames(bg.CCRF_CEM)[2] <- 'Score'
# bg.CCRF_CEM$Group='OtherSmallMolecules'
# bg.CCRF_CEM$Celline='CCRF_CEM'
# dim(bg.CCRF_CEM)
# head(bg.CCRF_CEM)
# 
# bg.MOLT_4<- bg.score[,c('Drug.name','LE.MOLT.4')]
# bg.MOLT_4 <- na.omit(bg.MOLT_4)
# bg.MOLT_4 <-   aggregate(bg.MOLT_4$LE.MOLT.4, by=list(type=bg.MOLT_4$Drug.name),mean)
# colnames(bg.MOLT_4)[2] <- 'Score'
# bg.MOLT_4$Group='OtherSmallMolecules'
# bg.MOLT_4$Celline='MOLT_4'
# dim(bg.MOLT_4)
# head(bg.MOLT_4)
# 

## ALL cell line
TALL.score <- NCI60.score[NCI60.score$Drug.name %in% candiDrug | 
                            NCI60.score$PubChem.SID %in% as.character(candiDrug.df$CID),
                          c('Drug.name','LE.CCRF.CEM','LE.MOLT.4')]
dim(TALL.score)
TALL.score$LE.CCRF.CEM <- as.numeric(TALL.score$LE.CCRF.CEM)
TALL.score$LE.MOLT.4 <- as.numeric(TALL.score$LE.MOLT.4)

TALL.CCRF_CEM<- TALL.score[,c('Drug.name','LE.CCRF.CEM')]
TALL.CCRF_CEM <- na.omit(TALL.CCRF_CEM)
TALL.CCRF_CEM <-   aggregate(TALL.CCRF_CEM$LE.CCRF.CEM, by=list(type=TALL.CCRF_CEM$Drug.name),mean)
colnames(TALL.CCRF_CEM)[2] <- 'Score'
TALL.CCRF_CEM$Group='CandidateDrugs'
TALL.CCRF_CEM$Celline='CCRF_CEM'
dim(TALL.CCRF_CEM)
head(TALL.CCRF_CEM)

TALL.MOLT_4<- TALL.score[,c('Drug.name','LE.MOLT.4')]
TALL.MOLT_4 <- na.omit(TALL.MOLT_4)
TALL.MOLT_4 <-   aggregate(TALL.MOLT_4$LE.MOLT.4, by=list(type=TALL.MOLT_4$Drug.name),mean)
colnames(TALL.MOLT_4)[2] <- 'Score'
TALL.MOLT_4$Group='CandidateDrugs'
TALL.MOLT_4$Celline='MOLT_4'
dim(TALL.MOLT_4)
head(TALL.MOLT_4)

### 最后靶向候选的药物（有NCI60的药物）
length(intersect(TALL.MOLT_4$type,TALL.CCRF_CEM$type))
TargetDrug <- intersect(TALL.MOLT_4$type,TALL.CCRF_CEM$type)
TargetDrug.df <- candiDrug.df[candiDrug.df$NAME %in% TargetDrug,]
write.table(TargetDrug.df,
            paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/targetDrug_',cutoff,'.txt',sep=''),
            col.names = T,sep='\t',row.names = F)


## rbind
# score.CCRF_CEM=rbind(TALL.CCRF_CEM,bg.CCRF_CEM)
# score.MOLT_4=rbind(TALL.MOLT_4,bg.MOLT_4)
# 
# score = rbind(score.CCRF_CEM,score.MOLT_4)

score = rbind(TALL.CCRF_CEM,TALL.MOLT_4)
score[1:2,]

plot.df <- score
colnames(plot.df)[1] <- 'Drug'
plot.df$DrugID <- candiDrugID[plot.df$Drug]
plot.df$Gene <- candiGene[plot.df$Drug]
plot.df$Tag <- paste(plot.df$DrugID,'\n(',plot.df$Drug,')',sep='')
head(plot.df)

#plot.df$Cellline=factor(plot.df$Cellline,levels = c('CCRF_CEM','MOLT_4'))

p.ALL <- ggplot() +  geom_bar(data = plot.df, 
                              aes(x = Tag, y = Score,fill=Celline),
                              stat = "identity",
                              width = 0.7, 
                              position = position_dodge(width = 0.9)) + 
  theme(panel.grid.major.x = element_line(colour = "black"), 
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size=8,colour = 'black'),
        legend.position = 'bottom') + 
  coord_flip()+
  scale_fill_manual(values = getPalette[2:1])+
  facet_wrap(~Gene,scales = "free_y")+
  labs(y="NCI60 Zscore")
p.ALL

ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/NCI60_',cutoff,'_v2.pdf',sep=''),
       p.ALL,width = 7,height = 5)
