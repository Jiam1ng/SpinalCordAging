library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stringr)
library(Hmisc)

all_gene = read.csv('E:/01_Project/04_SC/02_SC_bulk/DESeq2/SC_all_gene_anno.csv')
all_gene = all_gene[!duplicated(all_gene$Symbol),]

gene_list = all_gene$log2FoldChange
names(gene_list) = all_gene$Symbol
gene_list = sort(gene_list, decreasing = T)

go.bp = read.gmt('E:/04_Document/genelist/c5.go.bp.v2022.1.Hs.symbols.gmt')

res.bp = GSEA(gene_list, TERM2GENE = go.bp, eps = 0, pvalueCutoff = 0.05)
res.bp.df = res.bp@result
res.bp.df$pathway = str_sub(res.bp.df$Description, 6, nchar(res.bp.df$Description))
res.bp.df$pathway = tolower(res.bp.df$pathway)
res.bp.df$pathway = str_replace_all(res.bp.df$pathway, '_', ' ')
res.bp.df$pathway = capitalize(res.bp.df$pathway)
res.bp.df$logPadj = -log10(res.bp.df$p.adjust)

write.csv(res.bp.df, 'E:/01_Project/04_SC/02_SC_bulk/GSEA/GSEA_GO_BP.csv', row.names = F)

ggplot(res.bp.df, aes(NES, logPadj, size=logPadj)) +
  geom_point() +
  scale_x_continuous(breaks = c(-2,-1.5,-1,0,1,1.5,2)) +
  theme_classic() +
  theme(axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))

ggsave('E:/01_Project/04_SC/02_SC_bulk/GSEA/GSEA_GO_BP_point.pdf', height = 3.5, width = 5)

