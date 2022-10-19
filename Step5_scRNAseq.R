rm(list=ls())
library(Seurat)
library(ggplot2)
library(ggpubr)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)
show_col(getPalette)


plan='candiSet4'
cutoff='0'

###############
## 1.查看候选药物的靶点基因在TALL的哪些细胞中表达
###############
candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

pbmc.TALL <- readRDS('/local/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/tall_filtered_merged_umap.Rds')
DefaultAssay(pbmc.TALL) <- 'SCT'

Idents(pbmc.TALL) <- pbmc.TALL$class
CT <- levels(pbmc.TALL$class)
Color <- getPalette[1:length(CT)]

p0 <- VlnPlot(pbmc.TALL,features = candiGenes,group.by = 'class',pt.size = 0.001,
              ncol = ceiling(sqrt(length(candiGenes))),cols = Color)+
  theme(legend.position = 'none')
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/scRNAseqNC_candiGene_',cutoff,'_Vln.pdf',sep=''),
    width = 20,height = 20)
print(p0)
dev.off()



p1 <- DotPlot(pbmc.TALL,features = sort(candiGenes,decreasing = T),group.by = 'class')+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/scRNAseqNC_candiGene_',cutoff,'_Dot.pdf',sep=''),
    width = 13,height = 4)
print(p1)
dev.off()

### >25%和表达>1的基因
exp_info <- p1$data
exp_T <- exp_info[exp_info$id=='Malignant_T_Cells',]
exp_T$features.plot[exp_T$avg.exp>1&exp_T$pct.exp>25]


## 候选药物的靶点基因
candiGenes <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),sep='\t',header = T)
candiGenes <- unique(candiGenes$gene)
length(candiGenes)

p1 <- DotPlot(pbmc.TALL,features = candiGenes,group.by = 'class')+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))
p2 <- VlnPlot(pbmc.TALL,features = candiGenes,group.by = 'class',pt.size = 0.001,
              ncol = 5,cols = Color)+
  theme(legend.position = 'none')
# library(ggpubr)
# p <- ggarrange(p1,p2,widths = c(1, 1))

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/scRNAseqNC_candiDrugGene_',cutoff,'_Dot.pdf',sep=''),
       width = 8,height = 4)
print(p1)
dev.off()

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/scRNAseqNC_candiDrugGene_',cutoff,'_Vln.pdf',sep=''),
    width = 14,height = 7)
print(p2)
dev.off()

#### 只提取在Tumor T cell里表达超过>25% 的基因
candiSCGenes <- candiGenes[p1$data$pct.exp[1:length(candiGenes)]>25 &
             p1$data$avg.exp[1:length(candiGenes)]>1]
print(candiSCGenes)
saveRDS(candiSCGenes,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiSCGenes_',cutoff,'.RDS',sep=''))




###############
## 2. 只foucus到Tcell，查看Tcell发育阶段的marker
###############
rm(list=ls())
library(Seurat)
library(ggplot2)
library(scales)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/local/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/')
mtall <- readRDS('malignant_tall_subset.Rds')
# mtall$Cluster <- mtall$orig.ident
# mtall$Cluster[mtall$orig.ident=='SJTALL030263']='C1:ETP'
# mtall$Cluster[mtall$orig.ident=='SJTALL031201']='C2:ETP'
# mtall$Cluster[mtall$orig.ident=='X09_PATIENT']='C3:ETP'
# mtall$Cluster[mtall$orig.ident=='XB37']='C4:DP(CD3P)'
# mtall$Cluster[mtall$orig.ident=='XB41']='C5:SP (CD8)'
# mtall$Cluster[mtall$orig.ident=='XB47']='C6:SP (CD8)'
# mtall$Cluster <- factor(mtall$Cluster,
#                                levels =c('C1:ETP','C2:ETP','C3:ETP',
#                                          'C4:DP(CD3P)','C5:SP (CD8)','C6:SP (CD8)'))
# saveRDS(mtall,'malignant_tall_subset.Rds')

Color <- getPalette[1:length(unique(mtall$Cluster))]

p <- DimPlot(mtall,group.by = 'Cluster',cols = Color)
ggsave('umap_malignantT.pdf',p)


# VlnPlot(mtall,'CD44',split.by = 'Cluster')
# 
# VlnPlot(mtall,'CD38',split.by = 'Cluster')
# VlnPlot(mtall,'CD34',split.by = 'Cluster')
# VlnPlot(mtall,'CD1A',split.by = 'Cluster')
# 
# VlnPlot(mtall,'CD3D',split.by = 'Cluster')
# VlnPlot(mtall,'CD3E',split.by = 'Cluster')
# VlnPlot(mtall,'CD3G',split.by = 'Cluster')
# VlnPlot(mtall,'CD1E',split.by = 'Cluster')
# 
# VlnPlot(mtall,'CD4',split.by = 'Cluster')
# VlnPlot(mtall,'CCR7',split.by = 'Cluster')
# VlnPlot(mtall,'IL7R',split.by = 'Cluster')
# VlnPlot(mtall,'FAM102A',split.by = 'Cluster')
# VlnPlot(mtall,'CD28',split.by = 'Cluster')
# VlnPlot(mtall,'MAL',split.by = 'Cluster')
# VlnPlot(mtall,'ANK3',split.by = 'Cluster')
# VlnPlot(mtall,'AQP3',split.by = 'Cluster')
# 
# # CD8+ T
# VlnPlot(mtall,'CD8A',split.by = 'Cluster')
# VlnPlot(mtall,'CD8B',split.by = 'Cluster')
# VlnPlot(mtall,'PRDM1',split.by = 'Cluster')
# VlnPlot(mtall,'ZBTB38',split.by = 'Cluster')
# VlnPlot(mtall,'PRF1',split.by = 'Cluster')
# VlnPlot(mtall,'CST7',split.by = 'Cluster')
# 
# VlnPlot(mtall,'CD2',split.by = 'Cluster')
# VlnPlot(mtall,'CD7',split.by = 'Cluster')
# 
# VlnPlot(mtall,'CD19',split.by = 'Cluster')
# VlnPlot(mtall,'CD33',split.by = 'Cluster')
# 
# VlnPlot(mtall,'PTCRA',split.by = 'Cluster')
# VlnPlot(mtall,'CD27',split.by = 'Cluster')
# VlnPlot(mtall,'CD28',split.by = 'Cluster')
# 
# VlnPlot(mtall,'TRAC',split.by = 'Cluster')
# VlnPlot(mtall,'TRBC1',split.by = 'Cluster')
# VlnPlot(mtall,'TRDC',split.by = 'Cluster')
# VlnPlot(mtall,'TRGC1',split.by = 'Cluster')
# 
# VlnPlot(mtall,'IL2RA',split.by = 'Cluster')
# 
# VlnPlot(mtall,'LMO2',split.by = 'Cluster')
# VlnPlot(mtall,'LYL1',split.by = 'Cluster')
# VlnPlot(mtall,'HOXA9',split.by = 'Cluster')
# VlnPlot(mtall,'TLX3',split.by = 'Cluster')
# VlnPlot(mtall,'TLX1',split.by = 'Cluster')
# VlnPlot(mtall,'NKX2-1',split.by = 'Cluster')
# VlnPlot(mtall,'TAL1',split.by = 'Cluster')
# 
# FeaturePlot(mtall,'CD3D')
# FeaturePlot(mtall,c('LMO2','LYL1','HOXA9','TLX3','NKX2-1','TAL1'))

p_CD3D <- FeaturePlot(mtall,'CD3D',min.cutoff = 1.5,order = F)
p_CD3D
p_CD1E <- FeaturePlot(mtall,'CD1E',order = F)
p_CD1E
P_liuETP <- FeaturePlot(mtall,'ETP_Signature',min.cutoff = 0,order = F)
P_liuETP
P_zhangETP <- FeaturePlot(mtall,'ETPALL_Zhang_Signature',min.cutoff = 0,order = T)
P_zhangETP

P_liu.nonETP <- FeaturePlot(mtall,'Graux_Signature',min.cutoff = 0,order = T)
P_liu.nonETP
P_zhang.nonETP <- FeaturePlot(mtall,'TALL_Zhang_Signature',min.cutoff = 0,order = T)
P_zhang.nonETP
p <- CombinePlots(plots = list(p_CD3D,p_CD1E,P_liuETP,P_liu.nonETP,P_zhangETP,P_zhang.nonETP),ncol = 2)

ggsave('Tcell_MK_paper.pdf',p,height = 10,width =8)


p_MK <- FeaturePlot(mtall,c('CD44','CD34','CD4','CD8A','KIT','IL2RA','CD24','CD27','CD3D','CD1E'),ncol=2)
ggsave('Tcell_MK.pdf',p_MK,height = 20,width =10)

p_MK <- FeaturePlot(mtall,c('CD44','CD34','CD4','CD8A','KIT','IL2RA','CD24','CD27','CD3D','CD1E'),ncol=2)
ggsave('Tcell_MK.pdf',p_MK,height = 20,width =10)

p_MK <- FeaturePlot(mtall,c('CD4','CD8A'),ncol=2,order = T)
ggsave('Tcell_MK_1.pdf',p_MK,height = 4,width =10)


###############
## 3.查看候选药物的靶点基因在Tcell的表达阶段
###############
rm(list=ls())
library(Seurat)
library(ggplot2)
library(scales)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)

setwd('/local/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/')
mtall <- readRDS('malignant_tall_subset.Rds')
Color <- getPalette[1:length(unique(mtall$Cluster))]

plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

p0 <- VlnPlot(mtall,candiGenes,group.by = 'Cluster',pt.size = 0.0001,
              ncol = ceiling(sqrt(length(candiGenes))),cols=Color)+
  theme(legend.position = 'none')
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiGene_',cutoff,'_UMAP.pdf',sep=''),
    height = 20,width =20)
print(p0)
dev.off()


## 候选药物的靶点基因
candiGenes <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),sep='\t',header = T)
candiGenes <- unique(candiGenes$gene)
candiGenes <- candiGenes[c(2,6:9,1,3:5,10)]
length(candiGenes)

p1 <- VlnPlot(mtall,candiGenes,group.by = 'Cluster',pt.size = 0.0001,
              ncol = 5,cols=Color)
p2 <- FeaturePlot(mtall,candiGenes,ncol=2)
p3 <- DotPlot(mtall,features = candiGenes,group.by = 'Cluster')+
  theme(axis.text.x = element_text(angle = 30,hjust = 1))


pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiDrugGene_',cutoff,'_Vln.pdf',sep=''),
    height = 7,width =14)
print(p1)
dev.off()

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiDrugGene_',cutoff,'_FP.pdf',sep=''),
    height = 20,width =10)
print(p2)
dev.off()

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiDrugGene_',cutoff,'_Dot.pdf',sep=''),
    height = 4,width =8)
print(p3)
dev.off()
