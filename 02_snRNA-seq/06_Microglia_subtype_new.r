library(Seurat)
library(stringr)
library(dplyr)
library(ggpubr)
library(future)

sc.int = readRDS('SC_int_seu_obj.rds')
setwd('/data5/lijiaming/projects/01_single-cell/01_SC/07_Seurat_with_MN/04_Microglia/01_revise_final')
plan("multiprocess", workers = 5)

mic = subset(sc.int, celltype2 == 'Microglia')
mic = SCTransform(mic, verbose = T, variable.features.n = 500, assay='RNA')
mic = RunPCA(mic, verbose=F)

png('SC_mic_pca_elbow.png', height=500, width=500)
ElbowPlot(mic, ndims=50)
dev.off()

DefaultAssay(mic) = 'SCT'
mic = RunTSNE(mic, dims = c(1:20), seed.use=2022)
mic = FindNeighbors(mic, reduction = "pca", dims = c(1:20))
mic = FindClusters(mic, resolution=0.2, random.seed=2022)
mic$res0.2Ident = Idents(mic)

pdf('SC_mic_TSNE_cluster_res0.2.pdf', height = 4, width = 4)
TSNEPlot(mic, label=T, pt.size = 0.3) + NoLegend()
dev.off()

DefaultAssay(mic) = 'RNA'
mic = NormalizeData(mic)
cluster_markers = FindAllMarkers(mic, only.pos=T)
cluster_markers = subset(cluster_markers, p_val_adj<0.05 & avg_logFC>0.5)

write.csv(cluster_markers, 'SC_mic_cluster_markers_res0.2.csv', row.names=F)

for (gene in c('CSF1R', 'MRC1', 'CD163', 'CD247', 'MOBP', 'THEMIS','CHIT1', 'GPNMB')){
        p = FeaturePlot(mic, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('SC_mic_featureplot_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

png('SC_mic_featureplot_P2RY12.png', height = 500, width = 550, res = 100)
FeaturePlot(mic, features = 'P2RY12', raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
dev.off()

mic = NormalizeData(mic)
cluster.markers = FindAllMarkers(mic, only.pos=T)
cluster.markers = subset(cluster.markers, p_val_adj<0.05)

# exclude macrophage, oligodendrocyte and other non-microglia
mic.sub = subset(mic, res0.2Ident %in% c(0:2, 4, 5))
DefaultAssay(mic.sub) = 'RNA'
mic.sub = SCTransform(mic.sub, verbose = T, variable.features.n = 500)
mic.sub = RunPCA(mic.sub, verbose=F)

png('SC_mic_pca_subset_elbow.png', height=500, width=500)
ElbowPlot(mic.sub, ndims=50)
dev.off()

mic.sub = RunTSNE(mic.sub, dims = c(1:20), seed.use=2022)
mic.sub = FindNeighbors(mic.sub, reduction = "pca", dims = c(1:20))
mic.sub = FindClusters(mic.sub, resolution=0.2, random.seed=2022)
mic.sub$res0.2Ident = Idents(mic.sub)

pdf('SC_mic_subset_UMAP_cluster_res0.2.pdf', height = 4, width = 4)
DimPlot(mic.sub, label=T, pt.size=0.3) + NoLegend()
dev.off()

# exclude T cell
mic.sub = subset(mic.sub, res0.2Ident != '7')
DefaultAssay(mic.sub) = 'RNA'
mic.sub = SCTransform(mic.sub, verbose = T, variable.features.n = 500)
mic.sub = RunPCA(mic.sub, verbose=F)

png('SC_mic_pca_subset2_elbow.png', height=500, width=500)
ElbowPlot(mic.sub, ndims=50)
dev.off()

DefaultAssay(mic.sub) = 'SCT'
mic.sub = RunTSNE(mic.sub, dims = c(1:12), seed.use=2022)
mic.sub = FindNeighbors(mic.sub, reduction = "pca", dims = c(1:12))
mic.sub = FindClusters(mic.sub, resolution=0.2, random.seed=2022)
mic.sub$res0.2Ident = Idents(mic.sub)

pdf('SC_mic_subset2_UMAP_cluster_res0.2.pdf', height = 4, width = 4)
DimPlot(mic.sub, label=T, pt.size=0.3) + NoLegend()
dev.off()

DefaultAssay(mic.sub) = 'RNA'
mic.sub = NormalizeData(mic.sub)

for (gene in c('CSF1R', 'MRC1', 'CD163', 'CD247', 'MOBP', 'THEMIS','CHIT1', 'GPNMB')){
        p = FeaturePlot(mic.sub, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('SC_mic_subset2_featureplot_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

cluster_markers = FindAllMarkers(mic.sub, only.pos=T)
cluster_markers = subset(cluster_markers, p_val_adj<0.05 & avg_logFC>0.5)
write.csv(cluster_markers, 'SC_mic_subset2_cluster_markers_res0.2.csv', row.names=F)

dir.create('subset2_marker_feature')

top10 = cluster_markers %>% group_by(cluster) %>% top_n(n=10, wt=-p_val_adj)
for (i in 1:dim(top10)[1]){
        gene = as.character(top10$gene[i])
        p = FeaturePlot(mic.sub, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('subset2_marker_feature/SC_mic_subset2_featureplot_',i, '_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

# merge cluster 0, 1, and 3 for lack of specific marker genes
mic.sub = RenameIdents(mic.sub, c('0' = 'MG1', '1' = 'MG1', '3' = 'MG1',
                                  '2' = 'MG2', '4' = 'MG3', '5' = 'MG4'))
mic.sub$subtype = Idents(mic.sub)
subset_markers = FindAllMarkers(mic.sub, only.pos=T)
subset_markers = subset(subset_markers, p_val_adj<0.05 & avg_logFC>0.5)
write.csv(subset_markers, 'SC_mic_subset2_subtype_markers.csv', row.names=F)

dir.create('subtype_marker_feature')
top10 = subset_markers %>% group_by(cluster) %>% top_n(n=10, wt=-p_val_adj)
for (i in 1:dim(top10)[1]){
        gene = as.character(top10$gene[i])
        p = FeaturePlot(mic.sub, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('subtype_marker_feature/SC_mic_subset2_featureplot_',i, '_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

Idents(mic.sub) = mic.sub$subtype
pdf('SC_mic_subset2_TSNE_subtype.pdf', height = 4, width = 4)
TSNEPlot(mic.sub, label=T, pt.size = 0.3, cols=c('grey80', '#6bb9d2', '#d55640', '#469d88')) + NoLegend()
dev.off()

pdf('SC_mic_subset2_TSNE_subtype_noLabel.pdf', height = 4, width = 4)
TSNEPlot(mic.sub, pt.size = 0.3, cols=c('grey80', '#6bb9d2', '#d55640', '#469d88')) + NoLegend()
dev.off()

Idents(mic.sub) = mic.sub$age
pdf('SC_mic_subset2_TSNE_age_split.pdf', height = 4, width = 7)
TSNEPlot(mic.sub, label=F, pt.size = 0.3, split.by='age', cols = c('#b7cedb', '#48559a')) + NoLegend()
dev.off()


subset_markers$logPadj = -log10(subset_markers$p_val_adj)
subset_markers[subset_markers$logPadj>200,]$logPadj = 200
subset_markers$Rank = c(1:3, 1:25, 1:103, 1:7)
subset_markers = as.data.frame(subset_markers)
subset_markers$cluster = factor(subset_markers$cluster, levels=c('MG1', 'MG2', 'MG3', 'MG4'))

p = ggplot(subset_markers, aes(x=Rank, y=logPadj, color=cluster)) +
  geom_point(alpha=0.8) +
  facet_wrap(~cluster, nrow=1) +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_color_manual(values=c('##b7b7b7', '#6bb9d2', '#d55640', '#469d88')) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
        legend.position = 'none', axis.text = element_text(color='black'))

ggsave(plot=p, 'SC_mic_subset2_subtype_markers_point.pdf', height=4, width=6)

# cell proportion comparison
meta = mic.sub[[]]
cell.prop = data.frame(table(meta$subtype, meta$sample))
colnames(cell.prop) = c('subtype', 'sample', 'number')
cell.prop$age = factor(str_sub(cell.prop$sample, 1, 1), levels = c('Y', 'O'))
cell.prop = cell.prop %>% group_by(sample) %>% mutate(sampleSum = sum(number))
cell.prop$Cell.proportion = cell.prop$number / cell.prop$sampleSum * 100

p = ggplot(cell.prop, aes(age, Cell.proportion, fill = age)) +
      geom_boxplot(outlier.colour = NA, width = 0.6, color = 'black') +
      geom_jitter(width = 0.2) +
      scale_fill_manual(values = c('#b2cfdd', '#46549f')) +
      stat_compare_means(label='p.format', label.x=1.5, method='wilcox.test') +
      facet_wrap(~subtype, scales = 'free', nrow = 1) +
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
            legend.position = 'none', axis.text = element_text(color='black'))

ggsave(plot=p, filename='SC_mic_subset2_subtype_prop_boxplot.pdf', height = 3, width = 7)


total_prop = data.frame(table(meta$subtype))
colnames(total_prop) = c('subtype', 'number')
total_prop$group = 'total'
total_prop = total_prop[,c(1,3,2)]

age_prop = data.frame(table(meta$subtype, meta$age))
colnames(age_prop) = c('subtype', 'group', 'number')

total_prop = rbind(total_prop, age_prop)
total_prop = total_prop %>% group_by(group) %>% mutate(group_sum = sum(number))
total_prop$prop = total_prop$number / total_prop$group_sum * 100
total_prop$group = factor(total_prop$group, levels=c('total', 'Y', 'O'))
total_prop$subtype = factor(total_prop$subtype, levels=c('MG4', 'MG3', 'MG2', 'MG1'))

p = ggplot(total_prop, aes(group, prop, fill = subtype)) +
geom_bar(stat='identity', position='fill', width=0.8) +
scale_fill_manual(values=rev(c('#b7b7b7', '#6bb9d2', '#d55640', '#469d88'))) +
scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.1)) +
theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
            axis.text = element_text(color='black'), axis.text.x=element_text(angle=45, hjust=1, vjust=1))

ggsave(plot=p, filename='SC_mic_subset2_subtype_prop_bar.pdf', height = 8, width = 3)
