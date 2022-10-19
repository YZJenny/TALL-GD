rm(list=ls())
library("survival")
library("survminer")
library(survivalROC)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

## candiSet1: EMBO的DEG/ATAC都用作者EV9的
## candiSet2: EMBO的DEG用自己算的,ATAC用作者的EV9
## candiSet3: EMBO的DEG用自己算的,ATAC用作者的EV11
## candiSet4: EMBO的DEG用自己算的,ATAC自己算(EV7)

plan='candiSet4'
cutoff='0' 

### ATAC & EMBO
if(plan=='candiSet1'|plan=='candiSet2'){
  ATAC <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/GeneLst_paper.RDS')
  ATAC.up <- ATAC$up.Genes 
}else if(plan=='candiSet3'){
  EMBO_ATAC <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/GeneLst.RDS')
  ATAC.up <- EMBO_ATAC$up.Genes
}else if(plan=='candiSet4'){
  EMBO_ATAC <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/EMBO/GeneLst_EV7.RDS')
  ATAC.up <- na.omit(EMBO_ATAC$up.Genes)
}

print(length(ATAC.up))

### DEG
DEG0 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/EMBO/GeneLst.RDS')
DEG0.up <- DEG0$up.Genes
DEG0.dn <- DEG0$dn.Genes

DEG1 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE26713/GeneLst.RDS')
DEG1.up <- DEG1$up.Genes
DEG1.dn <- DEG1$dn.Genes

DEG2 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE48558/GeneLst.RDS')
DEG2.up <- DEG2$up.Genes
DEG2.dn <- DEG2$dn.Genes

DEG3 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DEG/GSE66638_GSE50999/GeneLst.RDS')
DEG3.up <- DEG3$up.Genes
DEG3.dn <- DEG3$dn.Genes

if(plan=='candiSet1'){
  ol.upDEG=intersect(intersect(DEG1.up,DEG2.up),ATAC.up)
}else if(plan=='candiSet2' | plan=='candiSet3' | plan=='candiSet4'){
  ol.upDEG=intersect(intersect(DEG1.up,DEG2.up),DEG0.up)
}
print(length(ol.upDEG))
print(length(intersect(ATAC.up,ol.upDEG)))

if(cutoff=='0'|cutoff=='0.5'){
  ### DepMap
  DepMap <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/GeneSet/DepMap/effectGene_',cutoff,'.RDS',sep=''))
  print(length(DepMap))
  
  candiGenes <- intersect(intersect(ATAC.up,ol.upDEG),DepMap)
  print(length(candiGenes))
  }else if(cutoff=='surv'){
  ### survival
  Surv <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/prognostic/TARGET/GeneLst.RDS')
  print(length(Surv))
  candiGenes <- intersect(intersect(ATAC.up,ol.upDEG),Surv)
  print(length(candiGenes))
  
  ## 根据surv图，判断哪些基因是riskgene
  ###################
  EXP <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/inputExp.txt',header = T,sep='\t')
  EFS <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/EFS_subtype.txt',header = T,row.names = 1,sep='\t')
  OS <- read.table('/local/yanzijun/CRU/TALL/data/NG_2017/OS_subtype.txt',header = T,row.names = 1,sep='\t')
  
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
    
    res.cut <- surv_cutpoint(data,time = "time", event = "status",variables = candiGenes)
    res.cat <- surv_categorize(res.cut)
    write.table(summary(res.cut),
                paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/Cut_',type,'.txt',sep=''),
                row.names = T,col.names = T,quote = F,sep='\t')
    
    write.table(res.cat,
                paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/Cat_',type,'.txt',sep=''),
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
      
      pvalue.plots[[k]] <- ggsurvplot(fit,data = clin_info, risk.table = TRUE,pval = T,
                                      palette = 'nejm',title = paste(type,' of ',Genes[k],sep=''))
      # pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/',Genes[k],'_',type,'_',cutoff,'.pdf',sep=''),
      #     width = 6,height = 6)
      # print(p, newpage = FALSE)
      # dev.off()
    }
    pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Surv/',type,'.pdf',sep=''),width = 6,height = 6)
    print(pvalue.plots)
    dev.off()
  }
  
  if(plan=='candiSet2'){
    protectGenes <- c('ATF3', 'TUBA1C')
  }else if(plan=='candiSet3'){
    protectGenes <- c('YWHAG','ZNF22','ATF3')
  }
  
  candiGenes <- candiGenes[!candiGenes %in% protectGenes]
}
print(length(candiGenes))
saveRDS(candiGenes,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))



###########################
### venn plot
###########################
## DEG
library(VennDiagram)
venn.diagram(x=list(Busra=DEG0.up,GSE26713=DEG1.up,GSE48558=DEG2.up), 
             filename=paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/venn_upDEG.png',sep=''), imagetype="png", 
             col="white", lwd=2,lty=2,
             fill=c(colors()[616], colors()[38], colors()[468]), alpha=c(0.6, 0.6, 0.6), 
             reverse=TRUE,
             height = 1600, width = 1600, resolution =300,  
             cat.cex=1.2,cex=1.5 #字号大小
)


## DAG, DEG & Dependent Genes
library(VennDiagram)

p=venn.diagram(x=list(Busra_ATAC=ATAC.up,
                    Busra_RNAseq=DEG0.up,GSE26713=DEG1.up,GSE48558=DEG2.up,
                    DepMap=DepMap), 
             filename=NULL,
             #filename=paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/venn_candiGenes.tiff',sep=''), imagetype="tiff", 
             col="white", lwd=2,lty=2,
             fill=getPalette[1:5], alpha=c(rep(0.5, 5)), 
             reverse=TRUE,
             height = 1800, width = 1800, resolution =300,  
             cat.cex=1,cex=1.2, #字号大小
             fontfamily = "serif",fontface = "bold",
             cat.fontface = "bold",cat.fontfamily = "serif"
             # cat.default.pos = "outer", # 位置, outer 内 text 外
             # cat.pos = c(0, 30, 150,-120,-30),  # 位置，用圆的度数
             # cat.dist = c(0.17, 0.002, 0.17, 0.17, 0.002),  # 位置，离圆的距离
)
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/venn_candiGenes.pdf',sep=''))
grid.draw(p)
dev.off()



### met_express
MET0 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/EMBO/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
MET1 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/GSE26713/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
MET2 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/GSE48558/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
MET3 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/GSE66638/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
#ol.MET=intersect(intersect(MET1,MET2),MET3)
ol.MET=intersect(MET1,MET2)
#ol.MET=intersect(MET0,MET1)
#ol.MET=MET0
print(length(ol.MET))


## met_express
venn.diagram(x=list(Busra=MET0,GSE26713=MET1,GSE48558=MET2), 
             filename="/local/yanzijun/CRU/TALL_FM/res/FIG/venn_enzyme.png", imagetype="png", 
             col="white", lwd=2,lty=2,
             fill=c(colors()[616], colors()[38], colors()[468]), alpha=c(0.6, 0.6, 0.6), 
             reverse=TRUE,
             height = 1600, width = 1600, resolution =300,  
             cat.cex=1.2,cex=1.5 #字号大小
             
)
