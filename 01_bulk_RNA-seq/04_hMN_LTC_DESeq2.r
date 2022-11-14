library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(DESeq2)

count.dir <- 'E:/00_Project/04_SC/08_MN_bulk_new/LTC/01_count/'

build_mat <- function(Count_dir){
  count_files <- dir(Count_dir)
  for (i in 1:length(count_files)){
    if (i==1){
      count_mat <- read.table(paste0(Count_dir, '/', count_files[i], '/', count_files[i], '_count.txt'), sep='\t', head=F)
      colnames(count_mat) <- c('Gene_ID', count_files[i])
    }else{
      tmp <- read.table(paste0(Count_dir, '/', count_files[i], '/', count_files[i], '_count.txt'), sep='\t', head=F)
      colnames(tmp) <- c('Gene_ID', count_files[i])
      count_mat <- merge(count_mat, tmp, by='Gene_ID')
    }
  }
  rownames(count_mat) <- count_mat$Gene_ID
  count_mat <- count_mat[6:dim(count_mat)[1],-1]
  return(as.matrix(count_mat))
}

samples = dir(count.dir)
all.raw <- build_mat(count.dir)

condition = c(rep('c1', 6), rep('c2', 5))
colData <- data.frame(samples, condition)

out_dir = 'E:/00_Project/04_SC/08_MN_bulk_new/LTC/02_DESeq2/'
prefix = 'SC_MN_LTC'
dir.create(out_dir)

geneID_symbol <- read.table('D:/Database/genome/hg19_gene_info_df.txt', sep='\t', header = F)
colnames(geneID_symbol) <- c('GeneID', 'Symbol', 'Gene_type')

dds <- DESeqDataSetFromMatrix(all.raw, colData, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
res <- results(dds2)
res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"GeneID"
resdata_anno<-merge(geneID_symbol, resdata, by="GeneID")
write.csv(resdata_anno, paste0(out_dir, prefix, "_all_gene_anno.csv"),row.names = F)

rld <- rlog(dds2)
write.csv(assay(rld), paste0(out_dir, prefix, "_rlog_data.csv"))

# PCA ####
rlog.exp = assay(rld)
pca.res = prcomp(t(rlog.exp))
summary(pca.res)
# PC1: 0.5476   , PC2: 0.1492

pca.data = as.data.frame(pca.res$x)
pca.data$condition = factor(c(rep('D18', 3), rep('D25', 3),
                              rep('D32', 3), rep('D40', 2)),
                            levels = c('D18', 'D25', 'D32', 'D40'))

loading = pca.res$rotation

pc1.gene = data.frame(PC1 = loading[,1])
pc1.gene$GeneID = rownames(loading)
pc1.gene = merge(geneID_symbol, pc1.gene, by ='GeneID')

write.csv(pc1.gene, paste0(out_dir, prefix, "_PC1_loadings.csv"),row.names = F)

pc1.gene = subset(pc1.gene, Gene_type == 'protein_coding')
pc1.gene = pc1.gene[order(pc1.gene$PC1, decreasing = T),]

pc1.gene.sub = subset(pc1.gene, abs(PC1) > 0.02 & abs(PC2) < 0.05)
write.csv(pc1.gene.sub, paste0(out_dir, prefix, "_PC1_DEG.csv"),row.names = F)
