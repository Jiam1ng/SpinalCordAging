library(Seurat)
library(stringr)
library(dplyr)


count.dir = '/data5/lijiaming/projects/01_single-cell/01_SC/MN_STRT/02_count/run/count/'
samples = dir(count.dir)
samples <- unlist(strsplit(dir(count.dir), "_count"))

samples <- samples[seq(1,length(samples),2)]

id.trans = read.table('/data5/lijiaming/data/genome/Macaca_fascicularis/MF5.0.94/MF_gene_info_v3.txt',
                      fill = T, sep='\t', header = T)
rownames(id.trans) = id.trans$GeneID
filtered.cell = read.csv('/data5/lijiaming/projects/01_single-cell/01_SC/MN_STRT/03_result/01_Seurat/MN_STRT_filtered_cell.csv')


ReadSeurat <- function(sample, qc.dir){
  tmp = read.table(paste0(count.dir, sample, '_count.csv'), 
                   sep = '\t', row.names = 1, header = T)
  tmp = tmp[rowSums(tmp) != 0,]
  tmp = as.matrix(tmp)
  rownames(tmp) = id.trans[rownames(tmp),]$Symbol
  
  tmp <- CreateSeuratObject(counts = tmp, min.cells = 3, min.features = 500, project = sample)
  tmp$batch <- sample
  tmp[["percent.mt"]] <- PercentageFeatureSet(tmp, pattern = "^MT-")
  
  tmp$cell = colnames(tmp)
  tmp$batch = sample
  
  p <- VlnPlot(tmp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  png(paste0(qc.dir, sample, '_qc_vln.png'), height = 500, width = 1500, res=100)
  print(p)
  dev.off()
  
  return(tmp)
}

sample_normalize <- function(tmp, sample, out.dir){
  plan("multiprocess", workers = 5)
  print(date())
  print(paste0(sample, ': SCTransform started'))
  tmp = subset(tmp, cell %in% as.character(filtered.cell$cell) & percent.mt < 15)
  tmp = SCTransform(tmp, verbose = FALSE)
  print(date())
  print(paste0(sample, ': Normalize finished'))
  
  return(tmp)
  print(paste0(sample ,' completed'))
}


for (sample in samples){
  tmp <- ReadSeurat(sample, '01_QC/MN_STRT_')
  assign(sample, tmp)
}

for (sample in samples){
  tmp <- sample_normalize(get(sample), sample, '01_QC/dims/MN_STRT_')
  assign(sample, tmp)
}

rm(list = setdiff(ls(), samples))
save.image('MN_STRT_samples_SCT.RData')

