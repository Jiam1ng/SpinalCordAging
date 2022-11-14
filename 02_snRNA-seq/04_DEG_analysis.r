library(dpylr)
library(ggplot2)
library(ggpubr)
library(stringr)

sc.int = readRDS('/data5/lijiaming/projects/01_single-cell/01_SC/04_Seurat_final/SC_int_seu_obj.rds')
celltypes = c("OPC", "OL", "Schwann", "Fib.Ast", "Pro.Ast", 
              "Microglia", "Ependymal", "Pericyte","EC", 
              'MN', "CNTN4.IN", "RBFOX1.IN")

deg = c()
Idents(sc.int) = paste(sc.int$celltype2, sc.int$age, sep = '_')
for (cell in celltypes)
  try({
    tmp = FindMarkers(sc.int, ident.1 = paste0(cell, '_O'), ident.2 = paste0(cell, '_Y'))
    tmp$gene = rownames(tmp)
    tmp$celltype = cell
    deg = rbind(deg, tmp)
  })

deg.sub = subset(deg, p_val < 0.05)

write.csv('SC_snRNAseq_DEG.csv', row.names=F)