library(Matrix)
library(future)
library(Seurat)
library(hdf5r)
library(stringr)
library(DoubletFinder)
library(dplyr)
library(biomaRt)

human <- useEnsembl("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", mirror="www")
mouse <- useEnsembl("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", mirror="www")

mm.trans <- getLDS(attributes = c('mgi_symbol', 'ensembl_gene_id'), mart = mouse,
                   attributesL = c('hgnc_symbol','ensembl_gene_id'), martL = human, uniqueRows=T)
colnames(mm.trans) <- c('mmGene', 'mmID','hsGene', 'hsID')

count = '/data5/lijiaming/projects/01_single-cell/01_SC/10_mouse/03_CellBender/v2/'
prefix = 'SC_mouse_'

samples <- c('OC-SC', 'YC-SC')

read10x_h5 = function (filename, use.names = TRUE, unique.features = TRUE) 
{
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = "r")
  genomes <- names(x = infile)
  output <- list()
  
  for (genome in genomes) {
    counts <- infile[[paste0(genome, "/data")]]
    indices <- infile[[paste0(genome, "/indices")]]
    indptr <- infile[[paste0(genome, "/indptr")]]
    shp <- infile[[paste0(genome, "/shape")]]
    #        features <- infile[[paste0(genome, "/", feature_slot)]][]
    
    features <- infile[[paste0(genome, "/features/name")]][]
    
    barcodes <- infile[[paste0(genome, "/barcodes")]]
    sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                               x = as.numeric(x = counts[]), dims = shp[], giveCsparse = FALSE)
    
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
    if (infile$exists(name = paste0(genome, "/features"))) {
      types <- infile[[paste0(genome, "/features/feature_type")]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(X = types.unique, FUN = function(x) {
          return(sparse.mat[which(x = types == x), ])
        }, simplify = FALSE, USE.NAMES = TRUE)
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  }
  else {
    return(output)
  }
}

dir.create('01_QC')
for (i in 1:length(samples)){
  h5 = paste0(count, samples[i], '_out_filtered.h5')
  tmp = read10x_h5(h5)
  
  tmp = tmp[intersect(rownames(tmp), mm.trans$mmGene),]
  rownames(tmp) = mm.trans[match(rownames(tmp), mm.trans$mmGene),]$hsGene
  tmp = tmp[!duplicated(rownames(tmp)) & rownames(tmp) != '',]

  tmp = CreateSeuratObject(counts = tmp, project = samples[i], min.cells = 3, min.features = 200)
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  tmp$sample = samples[i]
  tmp$age = factor(ifelse(str_detect(samples[i], 'Y'), 'Y', 'O'),
                   levels = c('Y', 'O'))
  
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  png(paste0('01_QC/', prefix, samples[i], '_qc_vln.png'), height = 500, width = 1500, res=100)
  print(p)
  dev.off()
  
  tmp <- RenameCells(tmp, add.cell.id = samples[i])
  
  assign(samples[i], tmp)
  print(paste0('--------', samples[i], '--------'))
}

sample_normalize <- function(tmp, sample){
  plan("multiprocess", workers = 5)
  print(date())
  print(paste0(sample, ': SCTransform started'))
  tmp <- subset(tmp, subset = percent.mt < 5 & nFeature_RNA > 200)
  tmp <- SCTransform(tmp, verbose = FALSE)
  print(date())
  print(paste0(sample, ': SCTransform finished'))
  
  tmp <- RunPCA(tmp, verbose=F)
  
  png(paste0('01_QC/dims/', prefix, sample, '_pca_heatmap.png'), height = 3000, width = 700)
  DimHeatmap(tmp, dims=1:30, cells=500, balanced=T)
  dev.off()
  
  p<- ElbowPlot(tmp, ndims = 30)
  png(paste0('01_QC/dims/', prefix, sample, '_ElbowPlot.png'), height = 600, width = 700)
  print(p)
  dev.off()
  
  return(tmp)
  print(paste0(sample ,' completed'))
}

dir.create('01_QC/dims/')
for (sample in samples){
  tmp <- sample_normalize(get(sample), sample)
  assign(sample, tmp)
}

for (sample in samples){
  print(sample)
  print(get(sample))
}

dims = list(c(1:30), c(1:30))
res = c(1, 1)

for (i in 1:length(samples)){
  tmp <- get(samples[i])
  tmp <- RunUMAP(tmp, dims = dims[[i]], verbose=F)
  tmp <- FindNeighbors(tmp, reduction = "pca", dims = dims[[i]])
  tmp <- FindClusters(tmp, res=res[i])
  tmp[["cluster"]] <- Idents(tmp)
  assign(samples[i], tmp)
  print(paste0(samples[i], ' completed'))
}

plan("multiprocess", workers = 1)
pK.df <- data.frame(matrix(nrow=0, ncol=2))
colnames(pK.df) <- c("Sample", "Optimal_pK")
detach('package:Seurat')
unloadNamespace('Seurat')
unloadNamespace('plotly')
for (i in 1:length(samples)){
  sweep.res.list <- paramSweep_v3(get(samples[i]), PCs = dims[[i]], sct = T, num.cores=5)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- arrange(bcmvn, desc(BCmetric))$pK[1]
  tmp <- data.frame(Sample=samples[i], Optimal_pK=pK)
  pK.df <- rbind(pK.df, tmp)
  print(bcmvn)
  print(paste0("--------------", samples[i], " completed (", i, "/2)--------------"))
}

write.csv(pK.df, 'SC_mouse_doubletfinder_pK.csv', row.names=F)

est.prop = c(0.076, 0.076)
doublet.prop = c()
for (i in 1:length(samples)) {
  pK.use <- as.numeric(as.character(pK.df$Optimal_pK[i]))
  tmp <- get(samples[i])
  homotypic.prop <- modelHomotypic(tmp@meta.data$cluster)
  nExp_poi <- round(est.prop[i]*length(colnames(tmp)))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  tmp <- doubletFinder_v3(tmp, PCs = dims[[i]], pN = 0.25, pK = pK.use, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  tmp[["doublet"]] <- tmp[[paste("DF.classifications_0.25", pK.use, nExp_poi.adj, sep="_")]]
  prop <- nExp_poi.adj/length(tmp@meta.data$cluster)
  prop.tmp <- data.frame(Sample=samples[i], Number=nExp_poi.adj, Doublet_prop=prop)
  doublet.prop <- rbind(doublet.prop, prop.tmp)
  assign(samples[i], tmp)
  print(paste0("--------------", samples[i], " completed (", i, "/2)--------------"))
}

plan("multiprocess", workers = 5)
for (i in samples){
  tmp <- get(i)
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- SCTransform(tmp, verbose = FALSE)
  assign(i, tmp)
  print(paste0(i, ' completed'))
}

rm(list = setdiff(ls(), c(samples, 'samples')))
save.image('SC_mouse_singlet_SCT.rds')

int.list = list()
for (tmp in samples){
  int.list = append(int.list, get(tmp))
}

names(int.list) = samples
options(future.globals.maxSize = 19000 * 1024^2)

int.features <- SelectIntegrationFeatures(object.list = int.list, nfeatures = 3000)

int.list <- lapply(X = int.list, FUN = function(x) {
  x <- RunPCA(x, features = int.features, verbose = FALSE)
})

int.list <- PrepSCTIntegration(object.list = int.list, anchor.features = int.features, verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = int.list, normalization.method = "SCT", 
                                      anchor.features = int.features, verbose = F,
                                      reduction = "rpca")

sc.mouse <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT", verbose = FALSE)

meta = sc.mouse[[]]
meta = meta[,1:8]
sc.mouse@meta.data = meta

DefaultAssay(sc.mouse) = 'integrated'
sc.mouse <- RunPCA(sc.mouse, verbose = FALSE)
dir.create('02_clustering')
png('02_clustering/SC_mouse_pca_heatmap.png', height=3000, width=1000)
DimHeatmap(sc.mouse, dims = 1:50, cells = 500, balanced = TRUE)
dev.off()

png('02_clustering/SC_mouse_pca_elbow.png', height=500, width=500)
ElbowPlot(sc.mouse, ndims=50)
dev.off()

#choose dimension as 1:25
sc.mouse <- RunTSNE(sc.mouse, dims = c(1:25), seed.use=20220210)

sc.mouse <- FindNeighbors(sc.mouse, reduction = "pca", dims = c(1:25))
sc.mouse = FindClusters(sc.mouse, resolution=2, random.seed=20220210)
sc.mouse$res2Ident = Idents(sc.mouse)

DefaultAssay(sc.mouse) = 'RNA'
sc.mouse = NormalizeData(sc.mouse)

Idents(sc.mouse) = sc.mouse$res2Ident
sc.mouse = RenameIdents(sc.mouse, c("0"="OL",
                                    "1"="OL",
                                    "2"="OL",
                                    "3"="OL",
                                    "4"="OL",
                                    "5"="Astrocyte",
                                    "6"="OL",
                                    "7"="OL",
                                    "8"="Microglia",
                                    "9"="Astrocyte",
                                    "10"="OL",
                                    "11"="Neuron",
                                    "12"="OL",
                                    "13"="OL",
                                    "14"="Schwann",
                                    "15"="OL",
                                    "16"="Meningeal",
                                    "17"="Schwann",
                                    "18"="Meningeal",
                                    "19"="EC",
                                    "20"="Astrocyte",
                                    "21"="Microglia",
                                    "22"="Microglia",
                                    "23"="OPC",
                                    "24"="Ependymal",
                                    "25"="Neuron",
                                    "26"="Pericyte",
                                    "27"="Astrocyte",
                                    "28"="Meningeal",
                                    "29"="Neuron",
                                    "30"="Microglia",
                                    "31"="Schwann",
                                    "32"="Meningeal",
                                    "33"="OL",
                                    "34"="Neuron",
                                    "35"="Neuron",
                                    "36"="Neuron",
                                    "37"="Neuron",
                                    "38"="Neuron",
                                    "39"="EC",
                                    "40"="Neuron",
                                    "41"="Neuron",
                                    "42"="Neuron",
                                    "43"="Schwann",
                                    "44"="OL",
                                    "45"="TTN+",
                                    "46"="Astrocyte",
                                    "47"="Neuron",
                                    "48"="Meningeal"))
celltypes = c("OPC", "OL", "Schwann", "Astrocyte", 
              "Microglia", "Ependymal", "Meningeal", "EC", "Pericyte",
              "Neuron", "TTN+")

sc.mouse$celltype = factor(Idents(sc.mouse), levels = celltypes)
Idents(sc.mouse) = sc.mouse$celltype

saveRDS(sc.mouse, 'SC_mouse_seu_obj.rds')


