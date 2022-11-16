library(Seurat)
library(stringr)
library(dplyr)
library(ggpubr)
library(future)

plan("multiprocess", workers = 5)
immune = readRDS('immune_cell_seu_obj.rds')

mic = state(immune, IC_type == 'Microglia')
mic = SCTransform(mic, verbose = T, variable.features.n = 500, assay='RNA')
mic = RunPCA(mic, verbose=F)

png('SC_mic_pca_elbow.png', height=500, width=500)
ElbowPlot(mic, ndims=50)
dev.off()

DefaultAssay(mic) = 'SCT'
mic = RunTSNE(mic, dims = c(1:15), seed.use=202211)
mic = FindNeighbors(mic, reduction = "pca", dims = c(1:15))
mic = FindClusters(mic, resolution=0.15, random.seed=202211)
mic$res0.15Ident = Idents(mic)

pdf('SC_mic_TSNE_cluster_res0.15.pdf', height = 4, width = 4)
TSNEPlot(mic, label=T, pt.size = 0.3) + NoLegend()
dev.off()

DefaultAssay(mic) = 'RNA'
mic = NormalizeData(mic)
cluster_markers = FindAllMarkers(mic, only.pos=T)
cluster_markers = subset(cluster_markers, p_val_adj<0.05 & avg_logFC>0.5)
write.csv(cluster_markers, 'SC_mic_cluster_markers_res0.15.csv', row.names=F)

mic = RenameIdents(mic, c('0' = "Microglia1", '1' = 'Microglia1', '4' = 'Microglia1', '5' = 'Microglia1', '6' = 'Microglia1',
                          '2' = 'Microglia2', '3' = 'Microglia3'))
mic$microglia_state = Idents(mic)

Idents(mic) = mic$microglia_state
pdf('SC_mic_TSNE_state.pdf', height = 4, width = 4)
TSNEPlot(mic, label=T, pt.size = 0.3, cols=c('grey80', '#6bb9d2', '#d55640')) + NoLegend()
dev.off()

pdf('SC_mic_TSNE_state_noLabel.pdf', height = 4, width = 4)
TSNEPlot(mic, pt.size = 0.3, cols=c('grey80', '#6bb9d2', '#d55640')) + NoLegend()
dev.off()

Idents(mic) = mic$age
pdf('SC_mic_TSNE_age_split.pdf', height = 4, width = 7)
TSNEPlot(mic, label=F, pt.size = 0.3, split.by='age', cols = c('#b7cedb', '#48559a')) + NoLegend()
dev.off()

Idents(mic) = mic$microglia_state
state_markers = FindAllMarkers(mic, only.pos=T)
state_markers = subset(state_markers, p_val_adj<0.05 & avg_logFC>0.5)
write.csv(state_markers, 'SC_mic_state_markers.csv', row.names=F)

state_markers$logPadj = -log10(state_markers$p_val_adj)
state_markers[state_markers$logPadj>200,]$logPadj = 200
state_markers$Rank = c(1:3, 1:24, 1:91)

p = ggplot(state_markers, aes(x=Rank, y=logPadj, color=cluster)) +
  geom_point(alpha=0.8, shape=16, size=2) +
  facet_wrap(~cluster, nrow=1) +
  scale_x_continuous(breaks = c(0,50,100)) +
  scale_color_manual(values=c('#b7b7b7', '#6bb9d2', '#d55640')) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
        legend.position = 'none', axis.text = element_text(color='black'))

ggsave(plot=p, 'SC_mic_state_markers_point.pdf', height=4, width=5.5)

saveRDS(mic, 'SC_microglia_seu_obj.rds')

# cell proportion comparison
meta = mic[[]]
cell.prop = data.frame(table(meta$subtype, meta$sample))
colnames(cell.prop) = c('subtype', 'sample', 'number')
cell.prop$age = factor(str_sub(cell.prop$sample, 1, 1), levels = c('Y', 'O'))
cell.prop = cell.prop %>% group_by(sample) %>% mutate(sampleSum = sum(number))
cell.prop$Cell.proportion = cell.prop$number / cell.prop$sampleSum * 100

p = ggplot(cell.prop, aes(age, Cell.proportion, fill = age)) +
      geom_boxplot(outlier.colour = NA, width = 0.6, color = 'black') +
      geom_jitter(width = 0.2) +
      scale_fill_manual(values = c('#b2cfdd', '#46549f')) +
      stat_compare_means(label='p.format', label.x=1.5, method='wilcox.test') +
      facet_wrap(~subtype, scales = 'free', nrow = 1) +
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
            legend.position = 'none', axis.text = element_text(color='black'))

ggsave(plot=p, filename='SC_mic_state2_subtype_prop_boxplot.pdf', height = 3, width = 7)


total_prop = data.frame(table(meta$subtype))
colnames(total_prop) = c('subtype', 'number')
total_prop$group = 'total'
total_prop = total_prop[,c(1,3,2)]

age_prop = data.frame(table(meta$subtype, meta$age))
colnames(age_prop) = c('subtype', 'group', 'number')

total_prop = rbind(total_prop, age_prop)
total_prop = total_prop %>% group_by(group) %>% mutate(group_sum = sum(number))
total_prop$prop = total_prop$number / total_prop$group_sum * 100
total_prop$group = factor(total_prop$group, levels=c('total', 'Y', 'O'))
total_prop$subtype = factor(total_prop$subtype, levels=c('MG4', 'MG3', 'MG2', 'MG1'))

p = ggplot(total_prop, aes(group, prop, fill = subtype)) +
geom_bar(stat='identity', position='fill', width=0.8) +
scale_fill_manual(values=rev(c('#b7b7b7', '#6bb9d2', '#d55640', '#469d88'))) +
scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.1)) +
theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
            axis.text = element_text(color='black'), axis.text.x=element_text(angle=45, hjust=1, vjust=1))

ggsave(plot=p, filename='SC_mic_state2_subtype_prop_bar.pdf', height = 8, width = 3)

# cell proportion comparison
meta = mic[[]]
cell.prop = data.frame(table(meta$microglia_state, meta$sample))
colnames(cell.prop) = c('state', 'sample', 'number')
cell.prop$age = factor(str_sub(cell.prop$sample, 1, 1), levels = c('Y', 'O'))
cell.prop = cell.prop %>% group_by(sample) %>% mutate(sampleSum = sum(number))
cell.prop$Cell.proportion = cell.prop$number / cell.prop$sampleSum * 100

p = ggplot(cell.prop, aes(age, Cell.proportion, fill = age)) +
      geom_boxplot(outlier.colour = NA, width = 0.6, color = 'black') +
      geom_jitter(width = 0.2) +
      scale_fill_manual(values = c('#b2cfdd', '#46549f')) +
      stat_compare_means(label='p.format', label.x=1.5, method='wilcox.test') +
      facet_wrap(~state, scales = 'free', nrow = 1) +
      theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
            legend.position = 'none', axis.text = element_text(color='black'))

ggsave(plot=p, filename='SC_mic_state_prop_boxplot.pdf', height = 3, width = 5)

total_prop = data.frame(table(meta$microglia_state))
colnames(total_prop) = c('state', 'number')
total_prop$group = 'total'
total_prop = total_prop[,c(1,3,2)]

age_prop = data.frame(table(meta$microglia_state, meta$age))
colnames(age_prop) = c('state', 'group', 'number')

total_prop = rbind(total_prop, age_prop)
total_prop = total_prop %>% group_by(group) %>% mutate(group_sum = sum(number))
total_prop$prop = total_prop$number / total_prop$group_sum * 100
total_prop$group = factor(total_prop$group, levels=c('total', 'Y', 'O'))
total_prop$state = factor(total_prop$state, levels=c('Microglia3', 'Microglia2', 'Microglia1'))

p = ggplot(total_prop, aes(group, prop, fill = state)) +
geom_bar(stat='identity', position='fill', width=0.8) +
scale_fill_manual(values=rev(c('#b7b7b7', '#6bb9d2', '#d55640'))) +
scale_y_continuous(expand=c(0,0), breaks=seq(0,1,0.1)) +
theme(panel.grid = element_blank(), panel.background = element_blank(),
            axis.line = element_line(color='black'), axis.ticks = element_line(color='black'), 
            axis.text = element_text(color='black'), axis.text.x=element_text(angle=45, hjust=1, vjust=1))

ggsave(plot=p, filename='SC_mic_state_prop_bar.pdf', height = 8, width = 3)
