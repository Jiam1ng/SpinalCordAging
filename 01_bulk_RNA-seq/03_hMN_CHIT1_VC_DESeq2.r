library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(DESeq2)

count.dir <- 'E:/00_Project/04_SC/08_MN_bulk_new/CHIT1/01_count'

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

geneID_symbol <- read.table('D:/Database/genome/hg19_gene_info_df.txt', sep='\t', header = F)
colnames(geneID_symbol) <- c('GeneID', 'Symbol', 'Gene_type')

samples = dir(count.dir)
all.count <- build_mat(count.dir)

# CHIT1 vs Control ####
compare <- 'CHIT1'
sample <- c(samples[1:6])
condition <- factor(group[1:6], levels=c('Control', 'CHIT1'))

count_mat <- all.count[,sample]
colData <- data.frame(sample, condition)

dds <- DESeqDataSetFromMatrix(count_mat, colData, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
res <- results(dds2)

rld <- rlog(dds2)
write.csv(assay(rld), paste0(out_dir, prefix, "_", compare, '_RNAseq_rlog_data.csv'))

res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"GeneID"

resdata_anno<-merge(geneID_symbol,resdata,by="GeneID")
write.csv(resdata_anno, paste0(out_dir, prefix, "_", compare, '_all_gene_anno.csv'),row.names = F)
diff_gene <-subset(resdata_anno, padj < 0.05 & abs(log2FoldChange) >= 1.5 & Gene_type == 'protein_coding')
write.csv(diff_gene, paste0(out_dir, prefix, "_", compare, '_diff_gene_anno.csv'),row.names = F)

# CHIT1+VC vs CHIT1 ####
compare <- 'CHIT1+VC'
sample <- c(samples[1:3], samples[10:12])
condition <- factor(c(group[1:3], group[10:12]), levels=c('CHIT1', 'VC+CHIT1'))

compare <- 'CHIT1'
sample <- c(samples[1:6])
condition <- factor(group[1:6], levels=c('Control', 'CHIT1'))

count_mat <- all.count[,sample]
colData <- data.frame(sample, condition)

dds <- DESeqDataSetFromMatrix(count_mat, colData, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
res <- results(dds2)

rld <- rlog(dds2)
write.csv(assay(rld), paste0(out_dir, prefix, "_", compare, '_RNAseq_rlog_data.csv'))

res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"GeneID"

resdata_anno<-merge(geneID_symbol,resdata,by="GeneID")
write.csv(resdata_anno, paste0(out_dir, prefix, "_", compare, '_all_gene_anno.csv'),row.names = F)
diff_gene <-subset(resdata_anno, padj < 0.05 & abs(log2FoldChange) >= 1.5 & Gene_type == 'protein_coding')
write.csv(diff_gene, paste0(out_dir, prefix, "_", compare, '_diff_gene_anno.csv'),row.names = F)
