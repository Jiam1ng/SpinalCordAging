library(dpylr)
library(ggplot2)
library(ggpubr)
library(stringr)

sc.int = readRDS('SC_int_seu_obj.rds')
celltypes = c("OPC", "OL", "Schwann", "Fib.Ast", "Pro.Ast", 
              "Microglia", 'Macrophage', 'T_cells', 'ImmunOL', 
              "Ependymal", "Pericyte","EC", 
              'MN', "CNTN4.IN", "RBFOX1.IN")

deg = c()
Idents(sc.int) = paste(sc.int$celltype_IC, sc.int$age, sep = '_')
for (cell in celltypes)
  try({
    tmp = FindMarkers(sc.int, ident.1 = paste0(cell, '_O'), ident.2 = paste0(cell, '_Y'))
    tmp$gene = rownames(tmp)
    tmp$celltype = cell
    deg = rbind(deg, tmp)
  })

deg.sub = subset(deg, p_val_adj < 0.05)

write.csv(deg.sub, 'SC_snRNAseq_DEG_with_immune_cell.csv', row.names=F)
