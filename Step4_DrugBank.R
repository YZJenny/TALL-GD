rm(list=ls())

###########################
### candiGenes的candiDrug, 手动搜索DrugBank网站 or 服务器文件（少）
###########################
plan='candiSet4'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)

######
## 手动搜索
######
wp=paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/Target/',sep='')
files <- list.files(wp)
candiDB <- c()
for(f in files){
  print(f)
  db <- read.csv(paste(wp,f,sep=''),header = T,sep='\t')
  db$Gene=f
  candiDB <- rbind(candiDB,db)
  print(all(colnames(db)==colnames(candiDB)))
}
candiDB$DETAILS <- NULL

### add InChI ID/Key
db2drg_out <- read.csv('/local/yanzijun/public/DrugBank/out.csv',header = T)
DBID2InChI_ID <- db2drg_out$InChI_ID
names(DBID2InChI_ID) <- db2drg_out$drugID

DBID2InChI_Key <- db2drg_out$InChI_Key
names(DBID2InChI_Key) <- db2drg_out$drugID

candiDB$InChI_ID <- DBID2InChI_ID[candiDB$DRUGBANK.ID]
candiDB$InChI_Key <- DBID2InChI_Key[candiDB$DRUGBANK.ID]
candiDB[1:3,]
write.table(candiDB,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),
            col.names = T,row.names = F,sep='\t',quote = F)



######
## 服务器文件
######  
# gene2db <- read.csv('/local/yanzijun/public/DrugBank/drug_target.txt',header = T,sep='\t')
# colnames(gene2db) <- c('gene','db')
# 
# db2drg <- read.csv('/local/yanzijun/public/DrugBank/drugbank_vocabulary.csv',header = T)
# print(length(unique(db2drg$DrugBank.ID))) #14594多种drug
# db2drg <- db2drg[,c(1,3)]
# colnames(db2drg) <- c('db','name')
# 
# db2drg_out <- read.csv('/local/yanzijun/public/DrugBank/out.csv',header = T)
# dim(db2drg_out);print(colnames(db2drg_out))
# db2drg_out <- db2drg_out[,c('drugID',"Drug_Category",'Drug_Type','Generic_Name','KEGGID',"CID","SID")]
# colnames(db2drg_out) <- c('db','category','type','name','keggid',"CID","SID")
# 
# drug.info <- merge(merge(db2drg,gene2db,all=TRUE),db2drg_out,all=TRUE)
# dim(drug.info)
# head(drug.info)
# 
# candiDB <- drug.info[drug.info$gene%in%candiGenes,]
# dim(candiDB)
# candiDrugs <- unique(candiDB$name)
# print(candiDrugs)
# print(table(candiDB$type))
# view(candiDB)
# write.table(candiDB,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiDrug_',cutoff,'.txt',sep=''),
#             col.names = T,row.names = F,sep='\t',quote = F)
# 
# ############
# ## candiDrug与GDSC中IC50的相关性
# ############
# .libPaths(c("/local/yzj/R/x86_64-pc-linux-gnu-library/4.0",
#             "/local/yanzijun/R/x86_64-pc-linux-gnu-library/4.0"))
# library(ggcorrplot)
# cor_data <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/GDSC/cor_',cutoff,'.RDS',sep=''))
# 
# new.df <- cor_data[rownames(cor_data) %in% candiDrugs,]
# print(intersect(rownames(cor_data),candiDrugs))
# if(length(intersect(rownames(cor_data),candiDrugs))>0){
#   new.df <- as.data.frame(new.df);colnames(new.df) <- intersect(rownames(cor_data),candiDrugs)
#   write.table(new.df,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/Drug/candiGDSC_',cutoff,'.txt',sep=''),
#               col.names = T,row.names = T,sep='\t',quote = F)
#   
#   p.GDSC <- pheatmap::pheatmap(t(new.df),cluster_cols = T,cluster_rows = F, 
#                                show_colnames = T,treeheight_col = 0,legend=FALSE,
#                                border_color = "white",
#                                colorRampPalette(c("navy", "white", "firebrick3"))(50),
#   )
#   ggsave(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/GDSC_',cutoff,'.pdf',sep=''),p.GDSC,height = 1.2,width = 7)
# }else{
#   print('no intersection!')
# }
# 
# 
# ###########################
# ### 从enzyme和C-map得到的结果
# ###########################
# ### candiDrug
# Drug0 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/cMap/EMBO/DrugLst.RDS')
# Drug1 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/cMap/GSE26713/DrugLst.RDS')
# Drug2 <- readRDS('/local/yanzijun/CRU/TALL_FM/res/GeneSet/cMap/GSE48558/DrugLst.RDS')
# candiDrugs <- unique(c(intersect(Drug0,Drug1),intersect(Drug0,Drug2),intersect(Drug1,Drug2),
#                        intersect(intersect(Drug0,Drug1),Drug2)))
# length(candiDrugs) #47
# print(intersect(intersect(Drug0,Drug1),Drug2))
# 
# ### candiMET
# MET0 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/EMBO/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
# MET1 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/GSE26713/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
# MET2 <- read.table('CRU/TALL_FM/res/GeneSet/met_express/GSE48558/enz_output_prediction_big0.txt',sep='\t',header = T)[,1]
# candiEnzyme=intersect(intersect(MET0,MET1),MET2)
# length(candiEnzyme) #76
# 
# candiDB <- db2drg$db[db2drg$name %in% c(intersect(db2drg$name,capitalize(candiDrugs)),
#   intersect(db2drg$name,candiDrugs))]
# 
# 
# 
# # print(length(unique(db2drg$cid)))
# # 
# cid2drg <- read.csv('/local/yanzijun/public/L1000/l10002cid.txt',header = F,sep='\t')
# colnames(cid2drg) <- c('name','cid')
# dim(cid2drg)
# candi_cid2drg <- cid2drg[cid2drg$name%in%candiDrugs,]
# # dim(candi_cid2drg)
# # db2drg <- merge(db2cid,cid2drg,by='cid')
# # head(db2drg)
# # candiDB <- db2drg$db[db2drg$cid %in% candi_cid2drg$cid]
# # length(candiDB)
# 
# 
# ## 看target
# candi.drugbank <- drugbank[drugbank$DrugBank_ID %in% candiDB,]
# dim(candi.drugbank)
# candi.drugbank <- candi.drugbank[candi.drugbank$X.GeneSymbol %in% candiEnzyme,]
# dim(candi.drugbank)
# 
# ## 看compound
# db2drg_out <- read.csv('/local/yanzijun/public/DrugBank/out.csv',header = T)
# dim(db2drg_out)
# db2drg_out <- db2drg_out[,c('drugID','CID','Generic_Name','KEGGID')]
# colnames(db2drg_out) <- c('db','cid','name','keggid')
# 
# candi.drugbank <- db2drg_out[db2drg_out$db %in% candiDB,]
# dim(candi.drugbank)
# 
# 
# load('/local/yanzijun/CRU/TALL_FM/scr/met_express/all_rn_gene_cpd.RData')
# gene2cp <- c()
# for(j in 1:length(candiEnzyme)){
#   enzyme=candiEnzyme[j]
#   rns <- gene_rn[[enzyme]]
#   compounds <- c()
#   for(i in 1:length(rns)){
#     rn=rns[i]
#     pro=all_rn_pro[[rn]]
#     sub=all_rn_sub[[rn]]
#     compounds <- c(compunds,c(pro,sub))
#   }
#   compounds <- unique(compounds)
#   #compounds <- paste0(compunds,';',collapse ='')
#   g2cp <- data.frame(gene=rep(enzyme,length(compounds)),compounds=compounds)
#   gene2cp <- rbind(gene2cp,g2cp)
# }
# 
# candiCP <- unique(candi.drugbank$keggid)
# candiCP <- candiCP[candiCP!='Not Available']
# candiCP
# 
# print(gene2cp[gene2cp$compounds %in% candiCP,])
# 
# 
# ## indirect PPI
# library(data.table)
# PPI <- fread('public/STRING/stringV10.5.txt',header = F)
# head(PPI)
# 
# length(candiEnzyme)
# 
# PPI_genes <- unique(c(PPI$V2[PPI$V1%in%candiEnzyme],PPI$V1[PPI$V2%in%candiEnzyme]))
# length(PPI_genes)
# 
# # 看target
# candi.drugbank <- drugbank[drugbank$DrugBank_ID %in% candiDB,]
# dim(candi.drugbank)
# candi.drugbank <- candi.drugbank[candi.drugbank$X.GeneSymbol %in% PPI_genes,]
# dim(candi.drugbank)
# length(unique(candi.drugbank$DrugBank_ID))

