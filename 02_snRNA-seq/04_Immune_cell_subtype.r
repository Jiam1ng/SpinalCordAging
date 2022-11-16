library(Seurat)
library(stringr)
library(dplyr)
library(ggpubr)
library(future)

sc.int = readRDS('SC_MN_int4-3_sub_seu.rds')

plan("multiprocess", workers = 5)

immune = subset(sc.int, celltype2 == 'Immune')
immune = SCTransform(immune, verbose = T, variable.features.n = 1000, assay='RNA')
immune = RunPCA(immune, verbose=F)

png('SC_immune_pca_elbow_1k.png', height=500, width=500)
ElbowPlot(immune, ndims=50)
dev.off()

DefaultAssay(immune) = 'SCT'
immune = RunTSNE(immune, dims = c(1:20), seed.use=2022)
immune = FindNeighbors(immune, reduction = "pca", dims = c(1:20))
immune = FindClusters(immune, resolution=0.2, random.seed=2022)
immune$res0.2Ident = Idents(immune)

pdf('SC_immune_1k_TSNE_cluster_res0.2.pdf', height = 4, width = 4)
TSNEPlot(immune, label=T, pt.size = 0.3) + NoLegend()
dev.off()

DefaultAssay(immune) = 'RNA'
immune = NormalizeData(immune)

for (gene in c('CSF1R', 'MRC1', 'CD163', 'CD247', 'MOBP', 'MOG', 'THEMIS', 'PTPRC')){
        p = FeaturePlot(immune, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('SC_immune_1k_featureplot_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

# exclude non-immune cells
immune = subset(immune, res0.2Ident != 10)
DefaultAssay(immune) = 'RNA'
immune = SCTransform(immune, verbose = T, variable.features.n = 800)
immune = RunPCA(immune, verbose=F)

png('SC_immune_pca_1k_immune_elbow.png', height=500, width=500)
ElbowPlot(immune, ndims=50)
dev.off()

immune = RunTSNE(immune, dims = c(1:22), seed.use=2022)
immune = FindNeighbors(immune, reduction = "pca", dims = c(1:22))
immune = FindClusters(immune, resolution=0.2, random.seed=2022)
immune$res0.2Ident = Idents(immune)

pdf('SC_immune_1k_immune_UMAP_cluster_res0.2.pdf', height = 4, width = 4)
DimPlot(immune, label=T, pt.size=0.3) + NoLegend()
dev.off()


DefaultAssay(immune) = 'RNA'
immune = NormalizeData(immune)

for (gene in c('CSF1R', 'MRC1', 'CD163', 'CD247', 'MOBP', 'MOG', 'THEMIS', 'PTPRC')){
        p = FeaturePlot(immune, features = gene, raster = F, pt.size=0.4, cols = c('grey90', '#244196'))
        png(paste0('SC_immune_1k_immune_featureplot_',gene,'.png'), height = 500, width = 550, res = 100)
        print(p)
        dev.off()
}

Idents(immune) = immune$res0.2Ident
immune = RenameIdents(immune, c('0' = 'Microglia', '1' = 'Microglia', '2' = 'Microglia',
                                '3' = 'Macrophage', '4' = 'Microglia', '5' = 'Microglia',
                                '6' = 'ImmunOL', '7' = 'Microglia', '8' = 'Microglia', '9' = 'T_cells'))

immune$IC_type = factor(Idents(immune), levels=c('Microglia', 'Macrophage', 'T_cells', 'ImmunOL'))

Idents(immune) = immune$IC_type
pdf('SC_immune_1k_immune_UMAP_celltype.pdf', height = 4, width = 4)
DimPlot(immune, label=T, pt.size=0.3, cols=c('#F78D3F', '#FCD271', '#317EDB', '#2BBBD8')) + NoLegend()
dev.off()

pdf('SC_immune_1k_immune_UMAP_celltype_noLabel.pdf', height = 4, width = 4)
DimPlot(immune, pt.size=0.3, cols=c('#F78D3F', '#FCD271', '#317EDB', '#2BBBD8')) + NoLegend()
dev.off()

markers = c('CSF1R', 'APBB1IP','MRC1', 'CD163', 'CD247', 'THEMIS', 'MOBP', 'MOG')
pdf('immune_cell_marker_dotplot.pdf', height = 3.5, width = 6.5)
DotPlot(immune, features = markers, cols = c('grey90', '#244196')) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
dev.off()

saveRDS(immune, 'immune_cell_seu_obj.rds')

immune.meta = immune[[]]
meta = sc.int[[]]

meta$celltype_IC = as.character(meta$celltype2)
meta[rownames(immune.meta),]$celltype_IC = as.character(immune.meta$IC_type)

celltypes = c("OPC", "OL", "Schwann", "Fib.Ast", "Pro.Ast", 
              "Microglia", 'Macrophage', 'T_cells', 'ImmunOL', 
              "Ependymal", 'Meningeal', "Pericyte","EC", 
              'MN', "CNTN4.IN", "RBFOX1.IN", "LGR5.IN")
meta$celltype_IC = factor(meta$celltype_IC, levels=celltypes)
sc.int@meta.data = meta

saveRDS(meta, 'SC_int_seu_with_IC_annot_meta.rds')

