rm(list=ls())
library(igraph)
library(parallel)
# 检测系统的CPU数
detectCores()

## ENSG to genename annotation file
# annoENSG <- read.table('/mdshare/node9/yzj/publicData/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gtf',header = F,sep='\t')
# ENSG2gene <- annoENSG$V9[annoENSG$V3=='gene']
# ENSG <- apply(as.matrix(ENSG2gene), 1, function(x) unlist(strsplit(unlist(strsplit(x,'; '))[1],'gene_id '))[2])
# gene_name <- apply(as.matrix(ENSG2gene), 1, function(x) unlist(strsplit(unlist(strsplit(x,'; '))[2],'gene_name '))[2])
# ENSG2gene <- data.frame(ENSG=ENSG,gene_name=gene_name)
# dim(ENSG2gene)
# write.table(ENSG2gene,'/mdshare/node9/yzj/publicData/annotation/hg19/Homo_sapiens.GRCh37.75.pc.ENSG2gene.txt',sep='\t',quote = F,row.names = F,col.names = T)

ENSG2gene <- read.table('/mdshare/node9/yzj/publicData/annotation/hg19/Homo_sapiens.GRCh37.75.pc.ENSG2gene.txt',sep='\t',header = T)

## candi Genes
plan='candiSet1'
cutoff='0'

candiGenes <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/candiGenes_',cutoff,'.RDS',sep=''))
length(candiGenes)
candiGenes.ENSG <- ENSG2gene$ENSG[ENSG2gene$gene_name %in% candiGenes]

## network
FunCoup <- read.table('/local/yanzijun/public/FunCoup5/FC5.0_H.sapiens_compact',header = F,sep='\t')
FunCoup <- FunCoup[,c(3,4,1)]
colnames(FunCoup) <- c('Gene1','Gene2','Weight')

## ALL-related genes
GAD <- read.csv('/local/yanzijun/public/GAD/GAD.txt',header = F,sep='\t')
TALL.GAD <- GAD$V1[grep('acute lymphoblastic leukemia|T-cell leukemia',GAD$V2,ignore.case = T)]
TALL.GAD.ENSG <- ENSG2gene$ENSG[ENSG2gene$gene_name %in% TALL.GAD]

g <- graph_from_data_frame(FunCoup,directed = FALSE)
path <- distances(g, v = V(g)[name %in% TALL.GAD.ENSG],to = V(g)[name %in% candiGenes.ENSG],
                  mode = 'all', weights =FunCoup$Weight,algorithm = "dijkstra")
mean_dist <- mean(path)

## permutation
func_perm <- function(x){
  print(x)
  random.ENSG <- sample(ENSG2gene$ENSG,length(candiGenes.ENSG))
  random.path <- distances(g, v = V(g)[name %in% TALL.GAD.ENSG],to = V(g)[name %in% random.ENSG],
                           mode = 'all', weights =FunCoup$Weight,algorithm = "dijkstra")
  mean_dist <- mean(random.path)
  return(mean_dist)
}

random_dist = unlist(mclapply(1:1000, func_perm, mc.cores = 40))

res <- list(mean_dist=mean_dist,random_dist=random_dist)
saveRDS(res,paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/GAD/dist_',cutoff,'.RDS',sep=''))


### plot
res <- readRDS(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/GAD/dist_',cutoff,'.RDS',sep=''))
mean_dist=res$mean_dist
random_dist=res$random_dist
mean_dist
pvalue <- length(which(random_dist<=mean_dist))/length(random_dist)
print(pvalue)

pdf(paste('/local/yanzijun/CRU/TALL_FM/res/',plan,'/FIG/dist_',cutoff,'.pdf',sep=''),4,4)
plot(density(random_dist),xlab='Distance',main='',lwd=2)
abline(v=mean_dist,col='red',lwd=2)
text(x=0.35,y=20,paste('P = ',pvalue,sep=''))
dev.off()
