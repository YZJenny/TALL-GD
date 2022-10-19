rm(list=ls())
library(Seurat)
library(ggplot2)
library(ggpubr)
library(scales)
library(future)
library(future.apply)

###############
## 1.MK
###############
plan("multisession", workers = 12) ###set the compute core
options(future.globals.maxSize = 60000 * 1024^2)

pbmc.TALL <- readRDS('/local/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/tall_filtered_merged_umap.Rds')
DefaultAssay(pbmc.TALL) <- 'SCT'

Idents(pbmc.TALL) <- pbmc.TALL$class
MK <- FindAllMarkers(pbmc.TALL,only.pos = T)
saveRDS(MK,"/local/yanzijun/CRU/scTALL/public/NC_2021/single_cell_rnaseq_input/markers.RDS")


