
setwd('E:/00_Project/04_SC/06_snRNA-seq_final/03_Seurat_MN/03_DEG/GO/')

celltypes = c("OPC", "OL", "Schwann", "Fib.Ast", "Pro.Ast", 
              "Microglia", "Ependymal", "Pericyte","EC", 
              'MN', "CNTN4.IN", "RBFOX1.IN")

deg = read.csv('E:/00_Project/04_SC/06_snRNA-seq_final/03_Seurat_MN/03_DEG/SC_snRNAseq_DEG.csv')
deg$DE = ifelse(deg$avg_logFC>0, 'up', 'down')

go.all = c()
for (cell in celltypes)
  try({
    go.raw = read.csv(paste0(cell, '.csv'))
    go.genes = as.character(go.raw$Symbols)
    go.genes = str_split(go.genes, '/')
    
    all.count = c()
    up.count = c()
    down.count = c()
    de.gene = c()
    de.count = c()
    for (gene.list in go.genes){
      tmp = length(gene.list)
      all.count = c(all.count, tmp)
      
      tmp = subset(deg, celltype == cell & DE == 'up' & gene %in% gene.list)
      up.gene = as.character(tmp$gene)
      up.count = c(up.count, length(up.gene))
      up.gene = paste(up.gene, collapse = '/')
      
      tmp = subset(deg, celltype == cell & DE == 'down' & gene %in% gene.list)
      down.gene = as.character(tmp$gene)
      down.count = c(down.count, length(down.gene))
      down.gene = paste(down.gene, collapse = '/')
      
      de.gene = rbind(de.gene, 
                      data.frame(Upregulated=up.gene,
                                 Downreguated=down.gene))
    }
    
    de.count = rbind(de.count, data.frame(upCount=up.count,
                                          downCount=down.count,
                                          allCount=all.count))
    
    go.out = cbind(go.raw, de.gene, de.count)
    colnames(go.out)[1] = 'GroupID'
    go.out$DE.score = (go.out$upCount - go.out$downCount) / (go.out$downCount + go.out$upCount)
    go.out$celltype = cell
    
    go.all = rbind(go.all, go.out)
    
    write.csv(go.out, paste0(cell, '_out.csv'), row.names = F)
  })

write.csv(go.all, 'E:/00_Project/04_SC/06_snRNA-seq_final/03_Seurat_MN/03_DEG/SC_GO_all_new.csv',row.names = F)
