library(RcisTarget)

deg = read.csv('SC_MN_CHIT1_diff_gene_anno.csv')
deg = deg[!str_detect(as.character(deg$Symbol), '^ENSMFAG'),]
up = subset(deg, log2FoldChange > 0)$Symbol
down = subset(deg, log2FoldChange < 0)$Symbol
geneLists = list(Up = as.character(up), Down = as.character(down))

data(motifAnnotations_hgnc)
motifRankings <- importRankings("/data2/zhengyandong/SCENIC_cisTarget_databases/hg19_cisTarget_databases/hg19-tss-centered-10kb-7species.mc9nr.feather")


motifs_AUC <- calcAUC(geneLists, motifRankings)

motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                           motifAnnot=motifAnnotations_hgnc)

motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                   geneSets=geneLists,
                                                   rankings=motifRankings, 
                                                   nCores=5,
                                                   method="iCisTarget")

split.df = c()
for (r in 1:dim(motifEnrichmentTable_wGenes)[1]){
  tmp = motifEnrichmentTable_wGenes[r,]
  tfs = tmp$TF_highConf
  tfs = unlist(strsplit(tfs, '; '))
  tmp = motifEnrichmentTable_wGenes[rep(r, length(tfs)),]
  tmp$TF.split = tfs
  split.df = rbind(split.df, tmp)
}

write.csv(split.df, 'SC_MN_CHIT1_RcisTarget_split.csv')
