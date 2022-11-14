library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(ggthemes)

count_dir <- 'E:/01_Project/04_SC/02_SC_bulk/count'
out_dir <- 'E:/01_Project/04_SC/02_SC_bulk/DESeq2'

geneID_symbol <- read.csv('D:/Database/01_Genome/02_Macaca_fascicularis_5.0/MF_gene_info_v2.txt', sep = '\t')

# build count matrix ####
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

sample <- factor(c("YF1", "YF2", "YF3", "YF4", "YM1", "YM2", "YM3", "YM4", 
                   "OF1", "OF2", "OF3", "OF4", "OM1", "OM2", "OM3", "OM4"), 
                 levels=c("YF1", "YF2", "YF3", "YF4", "YM1", "YM2", "YM3", "YM4", 
                          "OF1", "OF2", "OF3", "OF4", "OM1", "OM2", "OM3", "OM4"))
condition <- factor(c(rep('Young', 8), rep('Old', 8)), levels=c('Young', 'Old'))

count_mat <- build_mat(count_dir)
count_mat <- count_mat[,c(9:16,1:8)]
colnames(count_mat) <- sample

colData <- data.frame(sample, condition)

dds <- DESeqDataSetFromMatrix(count_mat, colData, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1,]
dds2 <- DESeq(dds)
res <- results(dds2)
rld <- rlog(dds2)
write.csv(assay(rld), paste0(out_dir,"/SC_RNAseq_rlog_data.csv"))

res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"GeneID"

resdata_anno<-merge(geneID_symbol,resdata,by="GeneID")
write.csv(resdata_anno, paste0(out_dir,"/SC_all_gene_anno.csv"),row.names = F)

diff_gene <-subset(resdata_anno, padj < 0.05 & abs(log2FoldChange) > 0.5)
write.csv(diff_gene, paste0(out_dir,"/SC_diff_gene_anno.csv"),row.names = F)

pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
pca_data$group <- condition
write.csv(pca_data,paste0(out_dir, 'SC_RNA-seq_PCA.csv'),row.names = F)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$gender = rep(c(rep('Female', 4), rep('Male', 4)), 2)
#4.5*3
#'grey80', '#AAD8B0', '#74A58E', '#725D2E'
ggplot(pca_data, aes(PC1, PC2, color=group, shape=gender)) +
  geom_point(size=3) +
  #  geom_text_repel(data=pca_data, aes(label=name),
  #                  color='black') +
  scale_color_manual(values=c('#b2cfdd', '#46549f')) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme_bw() +
  theme(axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))

ggsave(paste0(out_dir, "/SC_RNA-seq_PCA_gender.pdf"), height = 3, width = 4)


Volcano_data <- resdata_anno
Volcano_data$threshold <- factor(ifelse(Volcano_data$padj < 0.05 & abs(log2FoldChange) > 0.5,ifelse(Volcano_data$log2FoldChange > 0 ,'Upregulated','Downregulated'),'Unchanged'),
                                 levels = c('Upregulated', 'Unchanged', 'Downregulated'))

volcano_sum<-as.data.frame(t(summary(Volcano_data$threshold)))

Volcano_tile <- paste0('\nCutoff for pvalue is 0.05','\nThe number of upregulated gene is ',volcano_sum$`Upregulated`,'\nThe number of downregulated gene is ',volcano_sum$`Downregulated`,'\nThe number of unchanged gene is ',volcano_sum$Unchanged )


Volcano_data$lg10 = -log10(Volcano_data$padj)

#b21f1f
#1a2a6c

ggplot( data=Volcano_data,
        aes(x=log2FoldChange, y =lg10, colour=threshold, fill=threshold)) +
  scale_color_manual(values=c('#b89eb2', 'grey90', '#244196')) +
  geom_point(alpha=0.8, size=2.5) +
  xlim(c(-6.1,6.1)) +ylim(c(0,10))+
  ggtitle(Volcano_tile) +
  theme_few() +
  geom_vline(xintercept=c(-0.5,0.5),lty=2,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="grey",lwd=0.6)+
  #  theme(plot.title = element_text(hjust = 0.5,size=12),legend.title = element_blank())+
  labs(x="log2FC",y="-log10Padj") +
  theme(axis.text = element_text(color = 'black'))

ggsave(paste0(out_dir, "/SC_RNA-seq_DEG_volcano.pdf"), width = 6, height = 6.5)

