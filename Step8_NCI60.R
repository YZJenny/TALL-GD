rm(list=ls())

library(xlsx)
library(ggplot2)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)
show_col(getPalette)

## candi Genes
plan='candiSet4'
cutoff='0'

candiDrug.df <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),header = T,sep='\t')
candiDrug <- unique(candiDrug.df$NAME)

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

## bg: other cell line
bg.score <- NCI60.score[NCI60.score$Drug.name %in% candiDrug | 
                          NCI60.score$PubChem.SID %in% as.character(candiDrug.df$CID),
                        !colnames(NCI60.score) %in% c('PubChem.SID','LE.CCRF.CEM','LE.MOLT.4')]
dim(bg.score)
tmp=bg.score[,2:ncol(bg.score)]
tmp=as.data.frame(lapply(tmp,as.numeric))
tmp[is.na(tmp)] <- 0
bg.score=data.frame(Drug.name=bg.score$Drug.name,othercellline=apply(tmp,1,median))
#bg.score <- aggregate(othercellline~Drug.name, data=bg.score,FUN=mean)
length(unique(bg.score$Drug.name))

## ALL cell line
TALL.score <- NCI60.score[NCI60.score$Drug.name %in% candiDrug | 
                            NCI60.score$PubChem.SID %in% as.character(candiDrug.df$CID),
                          c('Drug.name','LE.CCRF.CEM','LE.MOLT.4')]
dim(TALL.score)
TALL.score$LE.CCRF.CEM <- as.numeric(TALL.score$LE.CCRF.CEM)
TALL.score$LE.MOLT.4 <- as.numeric(TALL.score$LE.MOLT.4)
#TALL.score <- aggregate(cbind(LE.CCRF.CEM,LE.MOLT.4)~Drug.name, data=TALL.score,FUN=mean)
TALL.score

#score=merge(TALL.score,bg.score,by='Drug.name')

score=cbind(TALL.score,bg.score)
score[,4] <- NULL
score[1:2,]
dim(TALL.score)
plot.df <- reshape2::melt(score)
colnames(plot.df)[1] <- 'Drug'
plot.df$DrugID <- candiDrugID[plot.df$Drug]
plot.df$Gene <- candiGene[plot.df$Drug]
plot.df$Tag <- paste(plot.df$DrugID,'\n(',plot.df$Drug,')',sep='')
head(plot.df)
colnames(plot.df) <- c('Drug','Cellline','Zscore','DrugID','Gene','Tag')
plot.df$Cellline=factor(plot.df$Cellline,levels = c('othercellline','LE.CCRF.CEM','LE.MOLT.4'))

p.ALL <- ggplot() +  geom_bar(data = plot.df, 
                     aes(x = Tag, y = Zscore,fill=Cellline),
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
  scale_fill_manual(values = c('grey',getPalette[2:1]))+
  facet_wrap(~Gene,scales = "free_y")+
  labs(y="NCI60 Zscore")
p.ALL
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/NCI60_',cutoff,'.pdf',sep=''),
       p.ALL,width = 7,height = 4)



p<- ggplot(plot.df, aes(x = Tag, y = Zscore, fill = Cellline))+
  geom_bar(stat = "summary", fun ="mean", position = position_dodge()) +
  stat_summary(fun.data = 'mean_sd', geom = "errorbar", colour = "black",
               width = 0.15,position = position_dodge( .8))+
  facet_wrap(~Gene,scales = "free_y")
q <- p+　geom_signif(comparisons = list(c('LE.CCRF.CEM','othercellline'),
                                         c('LE.MOLT.4','othercellline')))+
  theme_classic2()+
  coord_flip()+
  scale_fill_manual(values = c('grey',getPalette[2:1]))+
  facet_wrap(~Gene,scales = "free_y")+
  labs(x="")+
  theme(plot.title = element_text(hjust = 0.5,size = 10),
        legend.position = 'none')
q





## 所有的leukemia cell line
L.score <- NCI60.score[NCI60.score$Drug.name %in% candiDrug | NCI60.score$PubChem.SID %in% as.character(candiDrug.df$CID),
                          c(2,grep('^LE.',colnames(NCI60.score)))]
L.score$LE.CCRF.CEM <- as.numeric(L.score$LE.CCRF.CEM)
L.score$LE.MOLT.4 <- as.numeric(L.score$LE.MOLT.4)
L.score$LE.HL.60.TB. <- as.numeric(L.score$LE.HL.60.TB.)
L.score$LE.K.562 <- as.numeric(L.score$LE.K.562)
L.score$LE.RPMI.8226 <- as.numeric(L.score$LE.RPMI.8226)
L.score$LE.SR <- as.numeric(L.score$LE.SR)
L.score[is.na(L.score)] <- 0

L.score <- aggregate(cbind(LE.CCRF.CEM,LE.MOLT.4,LE.HL.60.TB.,LE.K.562,LE.RPMI.8226,LE.SR)~Drug.name, data=L.score,FUN=mean)
L.score

plot.df <- reshape2::melt(L.score)
colnames(plot.df) <- c('Drug','Cellline','Zscore')
plot.df$Cellline=factor(plot.df$Cellline)
p.L <- ggplot() +  geom_bar(data = plot.df, 
                          aes(x = Drug, y = Zscore,fill=Cellline),
                          stat = "identity",
                          width = 0.7, 
                          position = position_dodge(width = 0.9)) + 
  theme(panel.grid.major.x = element_line(colour = "black"), 
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank()) + 
  coord_flip()
p.L
ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/NCI60_',cutoff,'_L.pdf',sep=''),p.L,width = 6,height = 6)
