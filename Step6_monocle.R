rm(list=ls())
library(Seurat)
library(ggplot2)
library(monocle)
library(scales)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(9)
show_col(getPalette)

###############
## 查看候选单细胞药物的靶点基因在Tcell中的表达趋势
###############
plan='candiSet4'
cutoff='0'

# candiGenes <- read.table(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),sep='\t',header = T)
# candiGenes <- unique(candiGenes$gene)

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))

# candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiSCGenes_',cutoff,'.RDS',sep=''))
print(length(candiGenes))


min_expr=3
num_cells_expressed=30
qval=0.01
print(paste('min_expr=',min_expr,';num_cells_expressed=',num_cells_expressed,';qval=',qval,sep=''))
HSMM_root <- readRDS(paste('/local/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/monocle_SCT/monocle_',min_expr,'_',num_cells_expressed,'_',qval,'.RDS',sep=''))
cds_subset <- HSMM_root[candiGenes,]

CT <- pData(HSMM_root)$cluster

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiGene_',cutoff,'_SCT',min_expr,'_',num_cells_expressed,'_',qval,'_ps.pdf',sep=''),
    width = 13,height = 13)
# pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiSCGene_',cutoff,'_SCT',min_expr,'_',num_cells_expressed,'_',qval,'_ps.pdf',sep=''),
#     width = 13,height = 2.5)
plot_genes_in_pseudotime(cds_subset, color_by = "cluster",ncol = length(candiGenes)/6)+
    scale_fill_manual(values = getPalette[1:length(unique(pData(HSMM_root)$cluster))])+
    scale_color_manual(values = getPalette[1:length(unique(pData(HSMM_root)$cluster))])
dev.off()


## 每个分支单独画
upcells <- rownames(pData(cds_subset))[pData(cds_subset)$State %in% c('1','2')]
cds_subset_up <- cds_subset[,upcells]
table(pData(cds_subset_up)$State)
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiSCGene_',cutoff,'_SCT',min_expr,'_',num_cells_expressed,'_',qval,'_up_ps.pdf',sep=''),
    height = 2.5,width = 13)
plot_genes_in_pseudotime(cds_subset_up, color_by = "cluster",ncol = length(candiGenes))+
    scale_fill_manual(values = getPalette[1:length(unique(pData(HSMM_root)$cluster))])+
    scale_color_manual(values = getPalette[1:length(unique(pData(HSMM_root)$cluster))])
dev.off()

dncells <- rownames(pData(cds_subset))[pData(cds_subset)$State %in% c('1','3')]
cds_subset_dn <- cds_subset[,dncells]
table(pData(cds_subset_dn)$State)
pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiSCGene_',cutoff,'_SCT',min_expr,'_',num_cells_expressed,'_',qval,'_dn_ps.pdf',sep=''),
    height = 2.5,width = 13)
plot_genes_in_pseudotime(cds_subset_dn, color_by = "cluster",ncol = length(candiGenes))+
    scale_fill_manual(values = getPalette[1:length(unique(pData(HSMM_root)$cluster))])+
    scale_color_manual(values = getPalette[1:length(unique(pData(HSMM_root)$cluster))])
dev.off()


## 画基因的表达轨迹图
library(ggsci)
cds <- HSMM_root
colnames(pData(cds))
pData(cds)$CDK6 = log2( exprs(cds)['CDK6',]+1)
p1=plot_cell_trajectory(cds, color_by = 'CDK6')  + scale_color_gsea()
pData(cds)$PSMB2 = log2(exprs(cds)["PSMB2",]+1)
p2=plot_cell_trajectory(cds, color_by = "PSMB2") + scale_color_gsea()
pData(cds)$TUBA1A = log2(exprs(cds)['TUBA1A',]+1)
p3=plot_cell_trajectory(cds, color_by = 'TUBA1A') + scale_color_gsea()
pData(cds)$TUBB = log2(exprs(cds)['TUBB',]+1)
p4=plot_cell_trajectory(cds, color_by = 'TUBB') + scale_color_gsea()
pData(cds)$TYMS = log2(exprs(cds)['TYMS',]+1)
p5=plot_cell_trajectory(cds, color_by = 'TYMS') + scale_color_gsea()
library(patchwork)
p <- p1+p2+p3+p4+p5

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/candiSCGene_',cutoff,'_SCT',min_expr,'_',num_cells_expressed,'_',qval,'.pdf',sep=''),
    height = 6,width = 10)
print(p)
dev.off()
saveRDS(cds,paste('monocle_',type,'/monocle_root_',min_expr,'_',num_cells_expressed,'_',qval,'.RDS',sep=''))
