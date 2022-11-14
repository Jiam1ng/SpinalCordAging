library(Seurat)
library(stringr)
library(dplyr)
library(future)

load('/data5/lijiaming/projects/01_single-cell/01_SC/04_Seurat_final/SC_singlet_SCT_samples.RData')
load('/data5/lijiaming/projects/01_single-cell/01_SC/MN_STRT/03_result/01_Seurat/MN_STRT_samples_SCT.RData')


samples = c(setdiff(ls(), c("filtered.cell", "i", "meta", "tmp", "samples")))

int.list = list()
for (tmp in samples){
  int.list = append(int.list, get(tmp))
}

names(int.list) = samples
options(future.globals.maxSize = 19000 * 1024^2)

int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 2000)

int.list <- lapply(X = int.list, FUN = function(x) {
  x <- RunPCA(x, features = int.features, verbose = FALSE)
})

int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", 
                                      anchor.features = int.features, verbose = F,
                                      reduction = "rpca")

int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", 
                                      anchor.features = int.features, verbose = F,
                                      reduction = "rpca", k.anchor = 3)

sc.int <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

meta = sc.int[[]]
meta = meta[,c(colnames(meta)[1:8], 'batch', 'cell')]
rownames(filtered.cell) = filtered.cell$cell
meta[!is.na(meta$cell),]$sample = as.character(filtered.cell[as.character(meta[!is.na(meta$cell),]$cell),]$sample)
meta$age = factor(str_sub(meta$sample, 1, 1), levels = c('Y', 'O'))
meta$method = ifelse(is.na(meta$cell), '10x', 'STRT')

sc.int@meta.data = meta

sc.int <- RunPCA(sc.int, verbose = FALSE)
png('01_clustering/SC_pca_v4_heatmap.png', height=3000, width=1000)
DimHeatmap(sc.int, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png('01_clustering/SC_pca_v4_elbow.png', height=500, width=500)
ElbowPlot(sc.int, ndims=50)
dev.off()

DefaultAssay(sc.int) = 'integrated'
sc.int <- RunTSNE(sc.int, dims = c(1:30), seed.use=210817)
sc.int = FindClusters(sc.int, resolution=7, random.seed=210817)  
sc.int$res7Ident = Idents(sc.int)

DefaultAssay(sc.int) = 'RNA'
sc.int = NormalizeData(sc.int)

saveRDS(sc.int, 'SC_int_seu_obj.rds')
