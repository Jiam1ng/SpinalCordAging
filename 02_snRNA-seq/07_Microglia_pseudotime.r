library(Seurat)
library(monocle)

mic.sub = readRDS('/data5/lijiaming/projects/01_single-cell/01_SC/07_Seurat_with_MN/04_Microglia/SC_mic_subtype_seu.rds')
meta = mic.sub[[]]

exp_mat <- as.matrix(GetAssayData(mic.sub, slot='data'))
exp_mat <- exp_mat[rowSums(exp_mat)!=0,]
# cells = as.character(sample(colnames(exp_mat), 3000))
# sample_mat = exp_mat[, cells]
pd <- data.frame(mic.sub[[]])
fd <- data.frame(gene_short_name=rownames(exp_mat))
rownames(fd) <- fd$gene_short_name
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)

cds.obj <- newCellDataSet(as(exp_mat, "sparseMatrix"),
                          phenoData = pd,
                          featureData = fd)
cds.obj <- estimateSizeFactors(cds.obj)
cds.obj <- estimateDispersions(cds.obj)

marker.genes = read.csv('/data5/lijiaming/projects/01_single-cell/01_SC/07_Seurat_with_MN/04_Microglia/SC_mic_subtype_markers.csv')
marker.genes = subset(marker.genes, p_val_adj<0.05 & avg_logFC>0.5)
ordering_genes = unique(as.character(marker.genes$gene))

cds.obj <- setOrderingFilter(cds.obj, ordering_genes)
cds.obj <- reduceDimension(cds.obj.tmp, max_components = 2,
                           method = 'DDRTree')

cds.obj <- orderCells(cds.obj.tmp, root_state = '5')

saveRDS('monocle/SC_mic_monocle_obj.rds')

# branch DEG ####
BEAM_res <- BEAM(cds.obj, branch_point = 2, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

sig_gene <- subset(BEAM_res, qval < 1e-4)


plot_genes_branched_pseudotime(cds.obj[c('UCHL1', 'GRIA2'),],
                               branch_point = 2,
                               color_by = "subtype",
                               ncol = 1) +
  scale_color_manual(values=c('#f9beaf', '#4dbbd5', '#e64b35'))


branch_states = c('1', '2')
new_cds <- buildBranchCellDataSet(cds.obj[as.character(sig_gene$gene_short_name),], branch_states = branch_states,
        branch_point = 2, progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds.obj@dispFitInfo

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
    Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
    Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
BranchAB_exprs <- genSmoothCurves(new_cds[, ], cores = 15,
    trend_formula = "~sm.ns(Pseudotime, df=3) * Branch", relative_expr = T, new_data = rbind(newdataA,
        newdataB))  # time consuming

BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
    setdiff(pData(new_cds)$State, branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
    "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
    "Pseudotime"]))
BranchB_num <- BranchA_num

BranchA_exprs <- log10(BranchA_exprs + 1)
BranchB_exprs <- log10(BranchB_exprs + 1)

heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
    sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
    center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
    FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3

heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
    1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
exp_rng <- range(heatmap_matrix)
bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)

hmcols = colorRampPalette(c('#061a47','#97eaee', '#fccd3f', '#f04c46')) (length(bks) - 1)

hmcols = colorRampPalette(c('#46549f', '#c7e6f2', '#f9e1f6', '#b27da8')) (length(bks) - 1)

num_clusters = 5
ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE,
        cluster_rows = TRUE, show_rownames = F, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = "ward.D2",
        cutree_rows = num_clusters, silent = TRUE, filename = NA,
        breaks = bks, color = hmcols)

annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,
    num_clusters)))

branch_labels = c("Cell fate 1", "Cell fate 2")

colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
    `Cell Type` = c(rep(branch_labels[1], BranchA_num), rep("Pre-branch",
        2 * BranchP_num), rep(branch_labels[2], BranchB_num)))
colnames(annotation_col) <- "Cell Type"

branch_colors = c("#979797", "#F05662", "#7990C8")

names(branch_colors) <- c("Pre-branch", branch_labels[1],
        branch_labels[2])
annotation_colors = list(`Cell Type` = branch_colors)
names(annotation_colors$`Cell Type`) = c("Pre-branch", branch_labels)

feature_label <- row.names(heatmap_matrix)
row_ann_labels <- row.names(annotation_row)


row.names(heatmap_matrix) <- feature_label
row.names(annotation_row) <- row_ann_labels
ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE,
    cluster_rows = TRUE, show_rownames = F, show_colnames = F,
    clustering_distance_rows = row_dist, clustering_method = "ward.D2",
    cutree_rows = num_clusters, annotation_row = annotation_row,
    annotation_col = annotation_col, annotation_colors = annotation_colors,
    gaps_col = col_gap_ind, treeheight_row = 20, breaks = bks,
    fontsize = 6, color = hmcols, border_color = NA, silent = TRUE)

pdf('monocle/branch_DEG_hm.pdf', height = 6, width = 5)
print(ph_res)
dev.off()

deg.annot = annotation_row
deg.annot$gene = rownames(deg.annot)
colnames(sig_gene)[1] = 'gene'
deg.annot = merge(sig_gene, deg.annot, by = 'gene')
write.csv(deg.annot, 'monocle/branch_DEG_cluter.csv', row.names = F)

# branch expression trajectroy ####
library(reshape)
c2.genes = as.character(subset(deg.annot, Cluster == '2')$gene)
c2.exp = heatmap_matrix[c2.genes,]
c2.exp = melt(c2.exp)
colnames(c2.exp) = c('gene', 'cell', 'exp')
c2.exp$branch = ifelse(c2.exp$cell < 101, 'branch2', 'branch3')
c2.exp[c2.exp$cell > 100,]$cell = c2.exp[c2.exp$cell > 100,]$cell - 100
c2.exp$gene_group = paste(c2.exp$gene, c2.exp$branch, sep = '_')

p = ggplot(data = c2.exp, aes(x = cell,y = exp))  +
    geom_smooth(aes(group=gene_group, color = branch), se=F, method = 'loess') +
    scale_color_manual(values = c('#bde8ef', '#f7d4d0')) +
    geom_smooth(aes(group=branch), colour = "#D7281B", method = 'loess')+
    theme(panel.grid = element_blank(), panel.background = element_rect(color='black', fill='white'),
          axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = 'none')

pdf('monocle/cluter2_branch_fit.pdf', height = 4, width = 4)
print(p)
dev.off()


