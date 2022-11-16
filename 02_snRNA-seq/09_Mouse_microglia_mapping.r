library(Seurat)
library(stringr)
library(dplyr)
library(ggpubr)

mouse_sc = readRDS('SC_mouse_seu_obj.rds')
mouse_immune = subset(mouse_sc, celltype == 'Microglia')

DefaultAssay(mouse_immune) = 'RNA'
mouse_immune = SCTransform(mouse_immune, verbose = T, variable.features.n = 1000)
mouse_immune = RunPCA(mouse_immune, verbose=F)

png('mouse_immune_pca_elbow.png', height=500, width=500)
ElbowPlot(mouse_immune, ndims=50)
dev.off()

mouse_immune = RunTSNE(mouse_immune, dims = c(1:12), seed.use=2022)
mouse_immune = FindNeighbors(mouse_immune, reduction = "pca", dims = c(1:12))
mouse_immune = FindClusters(mouse_immune, resolution=0.2, random.seed=2022)
mouse_immune$res0.2Ident = Idents(mouse_immune)

pdf('mouse_immune_tSNE_cluster_res0.2.pdf', height = 4, width = 4)
TSNEPlot(mouse_immune, label=T, pt.size=1) + NoLegend()
dev.off()

DefaultAssay(mouse_immune) = 'RNA'
mouse_immune = NormalizeData(mouse_immune)

for (gene in c('CSF1R', 'MRC1', 'CD163', 'CD247', 'MOBP', 'THEMIS', 'PTPRC')){
        p = FeaturePlot(mouse_immune, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('mouse_mic_immune_featureplot_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

# microglia subset
mouse_mic = subset(mouse_immune, res0.2Ident %in% c(0, 1, 2))
DefaultAssay(mouse_mic) = 'RNA'
mouse_mic = SCTransform(mouse_mic, verbose = T, variable.features.n = 500)
mouse_mic = RunPCA(mouse_mic, verbose=F)

png('mouse_mic_pca_elbow.png', height=500, width=500)
ElbowPlot(mouse_mic, ndims=50)
dev.off()

mouse_mic = RunTSNE(mouse_mic, dims = c(1:12), seed.use=2022)
mouse_mic = FindNeighbors(mouse_mic, reduction = "pca", dims = c(1:12))
mouse_mic = FindClusters(mouse_mic, resolution=0.2, random.seed=2022)

# monkey microglia mapping
monkey_mic = readRDS('SC_microglia_seu_obj.rds')
DefaultAssay(monkey_mic) = 'SCT'
mic_anchors = FindTransferAnchors(reference = monkey_mic, query = mouse_mic, dims = 1:12)
mapping_res = TransferData(anchorset = mic_anchors, refdata = monkey_mic$microglia_state, dims = 1:12)
mouse_mic = AddMetaData(mouse_mic, metadata = mapping_res)

Idents(mouse_mic) = mouse_mic$predicted.id
pdf('mouse_mic_TSNE_state.pdf', height = 4, width = 4)
TSNEPlot(mouse_mic, label=T, pt.size = 1, cols=c('#d55640', 'grey80')) + NoLegend()
dev.off()

pdf('mouse_mic_TSNE_state_noLabel.pdf', height = 4, width = 4)
TSNEPlot(mouse_mic, pt.size = 1, cols=c('#d55640', 'grey80')) + NoLegend()
dev.off()

Idents(mouse_mic) = mouse_mic$age
pdf('mouse_mic_TSNE_age_split.pdf', height = 4, width = 7)
TSNEPlot(mouse_mic, label=F, pt.size = 1, split.by='age', cols = c('#48559a', '#b7cedb')) + NoLegend()
dev.off()

DefaultAssay(mouse_mic) = 'RNA'
mouse_mic = NormalizeData(mouse_mic)
mouse_markers = FindAllMarkers(mouse_mic, only.pos=T)
mouse_markers = subset(mouse_markers, p_val_adj<0.05 & avg_logFC>0.5)
write.csv(mouse_markers, 'mouse_state_markers.csv', row.names=F)

saveRDS(mouse_mic, 'mouse_mic_seu_obj.rds')

