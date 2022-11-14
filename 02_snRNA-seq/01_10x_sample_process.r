library(Matrix)
library(future)
library(Seurat)
library(hdf5r)
library(stringr)
library(DoubletFinder)
library(dplyr)

count = '/data5/lijiaming/projects/01_single-cell/01_SC/01_count/cellbender/'
qc.dir = '01_QC/'
prefix = 'SC_final_'

samples <- unlist(str_split(dir(count), "_"))
samples <- samples[seq(1,length(samples),3)]

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


for (i in 1:length(samples)){
  h5 = paste0(count, samples[i], '_out_filtered.h5')
  tmp = read10x_h5(h5)
  tmp = CreateSeuratObject(counts = tmp, project = samples[i], min.cells = 3, min.features = 200)
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  tmp$sample = samples[i]
  tmp$age = factor(ifelse(str_detect(samples[i], 'Y'), 'Y', 'O'),
                   levels = c('Y', 'O'))
  
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  png(paste0(qc.dir, prefix, samples[i], '_qc_vln.png'), height = 500, width = 1500, res=100)
  print(p)
  dev.off()
  
  tmp <- RenameCells(tmp, add.cell.id = samples[i])
  
  assign(samples[i], tmp)
  print(paste0('--------', samples[i], '--------'))
}

sample_normalize <- function(tmp, sample, out.dir){
  plan("multiprocess", workers = 5)
  print(date())
  print(paste0(sample, ': SCTransform started'))
  tmp <- subset(tmp, subset = percent.mt < 5 & nFeature_RNA > 200)
  tmp <- SCTransform(tmp, verbose = FALSE)
  print(date())
  print(paste0(sample, ': SCTransform finished'))
  
  tmp <- RunPCA(tmp, verbose=F)
  
  png(paste0(out.dir, prefix, sample, '_pca_heatmap.png'), height = 3000, width = 700)
  DimHeatmap(tmp, dims=1:30, cells=500, balanced=T)
  dev.off()
  
  p<- ElbowPlot(tmp, ndims = 30)
  pdf(paste0(out.dir, prefix, sample, '_ElbowPlot.pdf'), height = 6, width = 7)
  print(p)
  dev.off()
  
  return(tmp)
  print(paste0(sample ,' completed'))
}

for (sample in samples){
  tmp <- sample_normalize(get(sample), sample, 
                          '01_QC/dims/')
  assign(sample, tmp)
}


dims = list(c(1:12), c(1:10), c(1:8), c(1:7), 
            c(1:6), c(1:9), c(1:8), c(1:6), 
            c(1:7), c(1:11), c(1:7), c(1:7), 
            c(1:11), c(1:10), c(1:11), c(1:10))
res = c(0.7, 0.6, 0.5, 0.8, 0.6, 0.6, 0.6, 0.7, 
        0.7, 0.8, 0.9, 0.6, 0.9, 0.4, 0.4, 0.6)

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
  print(paste0("--------------", samples[i], " completed (", i, "/16)--------------"))
}

est.prop = c(0.054, 0.046, 0.039, 0.061, 
             0.05, 0.05, 0.054, 0.061,
             0.054, 0.061, 0.069, 0.061,
             0.069, 0.033, 0.031, 0.046)
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
  print(paste0("--------------", samples[i], " completed (", i, "/16)--------------"))
}

for (i in samples){
  tmp <- DimPlot(get(i), group.by='doublet', cols=c('firebrick3', 'grey90'))
  png(paste0('01_QC/doublet/SC_', i, '_UMAP_doublet.png'), height = 400, width = 450)
  print(tmp)
  dev.off()
}

rm(list = setdiff(ls(), samples))
save.image('SC_doublet_annot_samples.RData')

plan("multiprocess", workers = 5)
for (i in samples){
  tmp <- get(i)
  tmp <- subset(tmp, doublet=='Singlet')
  tmp <- SCTransform(tmp, verbose = FALSE)
  assign(i, tmp)
  print(paste0(i, ' completed'))
}

rm(list=setdiff(ls(), c(samples, 'samples')))
save.image('SC_singlet_SCT_samples.RData')

