library(DEXSeq)
library(DESeq2)
library(goseq)
library(limma)
library(genefilter)
library(e1071)

library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GO.db)
library(KEGG.db)
library(pathview)

library(tidyverse)
library(reshape2)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(DiagrammeR)
library(stringr)
library(knitr)

load('processed_data.RData')

source('R/00.functions.R')

# figure_1B
assay(dds) %>%
  melt %>%
  setNames(c('gene_id', 'sample', 'count')) %>%
  left_join(sample_table) %>%
  filter(gene_id %in% c('Gapdh', 'Lpl', 'Cebpa', 'Pparg'),
         stage %in% c('Proliferation', 'Differentiation')) %>%
  group_by(gene_id, stage) %>%
  summarise(ave = log(mean(count))) %>%
  ungroup %>%
  mutate(gene_id = factor(gene_id, levels = c('Gapdh', 'Lpl', 'Cebpa', 'Pparg'))) %>%
  ggplot(aes(x = gene_id, y = ave, fill = stage)) +
  geom_col(position = 'dodge', width = .5) +
  theme_classic() +
  labs(x = '', y = 'Average Log Count') +
  theme(legend.position = c(.8, 1),
        legend.direction = 'vertical', 
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(.5, units = 'cm'),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, size = 10, hjust = 1, face = 'italic'),
        plot.margin = unit(c(1,.5,.8,.2), 'cm'))
ggsave('manuscript/figures/figure1B.png', width = 4, height = 4, units = 'in', dpi = 500)

# figure_1C
ind <- c('Fasn', 'Acly', 'Acaca', 'Elovl6', 'Scd1', 'Scd2', 'Scd3', 'Scd4', 'Dgat1', 'Dgat2', 'Gapdh')
lab <- c('Scd4', 'Elovl6', 'Gapdh','Scd3', 'Acaca', 'Dgat1', 'Acly', 'Dgat2', 'Fasn', 'Scd2', 'Scd1')
assay(dds) %>%
  melt %>%
  setNames(c('gene_id', 'sample', 'count')) %>%
  left_join(sample_table) %>%
  filter(gene_id %in% ind) %>%
  group_by(gene_id, stage) %>%
  summarise(ave = log(mean(count))) %>%
  ggplot(aes(x = stage, y = ave, group = gene_id, color = gene_id)) +
  geom_line() +
  theme_classic() +
  labs(x = '', y = 'Average Log Count') +
  theme(legend.position = 'none',
        legend.text = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.margin = unit(c(1,.5,0,.2), 'cm')) +
  annotate('text', x = 4.1, y = c(7.2, 8.2, 9, 9.8, 10.6, 11.2, 11.8, 12.4, 13.4, 14, 15.8), label = lab, size = 3, hjust = 0, fontface = 'italic')

ggsave('manuscript/figures/figure1C.png', width = 4, height = 4, units = 'in', dpi = 500)

# figure_2a
png('manuscript/figures/figure2A.png', width = 3, height = 3, units = 'in', res = 500)
ind <- gene_annotation %>%
  filter(GO == 'GO:0006914') %>%
  dplyr::select(SYMBOL) %>%
  unlist %>%
  unique

res_groups %>%
  filter(row %in% ind) %>%
  mutate(stage = str_sub(stage, start = 6),
         padj = ifelse(padj < .1, 1, 0)) %>%
  acast(row ~ stage, value.var = 'padj') %>%
  modified_vennDiagram(circle.col = c('red', 'blue', 'green'),
                       mar = rep(0, 4),
                       cex = rep(.65, 4),
                       names = c('Clonal Expansion', 'Differentiation', 'Quiescence'))
dev.off()

png('manuscript/figures/figure2B.png', width = 3, height = 3, units = 'in', res = 500)
res_exon_groups %>%
  bind_rows %>%
  mutate(gene_id = str_sub(strsplit2(groupID, '_')[,2], end = -2),
         id = paste(groupID, featureID),
         padj = ifelse(padj < .1, 1, 0)) %>%
  filter(gene_id %in% ind) %>%
  acast(id ~ stage, value.var = 'padj') %>%
  modified_vennDiagram(circle.col = c('red', 'blue', 'green'),
                       mar = rep(0, 4),
                       cex = rep(.65, 4),
                       names = c('Clonal Expansion', 'Differentiation', 'Quiescence'))
dev.off()

# figure_2
# get autophagy gene set symbols
ind <- gene_annotation %>%
  filter(GO == 'GO:0006914') %>%
  dplyr::select(SYMBOL) %>%
  left_join(res_groups, by = c('SYMBOL'='row')) %>%
  filter(padj < .1, abs(log2FoldChange) > 1) %>%
  dplyr::select(SYMBOL) %>%
  unlist %>% unique

## figure: heatmap of significant autophagy genes
mat <- assay(dds)[ind,]
cols <- colorRampPalette(brewer.pal(n = 5, name = "Greens"))(200)
ann <- data.frame(Induction = sample_table$induction,
                  Stage = sample_table$stage,
                  row.names = sample_table$sample)
ann_cols <- colorRampPalette(brewer.pal(n = 4, name = 'Blues'))(4) 
names(ann_cols) <- unique(ann$`Time Point`)
ann_cols <- list(Induction = c(Non = '#FEE6CE', MDI = '#E6550D'),
                 Stage = c(Proliferation = '#EFF3FF',
                           Quiescence = '#BDD7E7',
                           `Clonal expansion` = '#6BAED6',
                           Differentiation = '#2171B5'))

png('manuscript/figures/figure2A.png', width = 4.5, height = 7, units = 'in', res = 500)
pheatmap(mat,
         scale = 'row', legend = FALSE,
         cellheight = 12, cellwidth = 15, fontsize = 8, cluster_cols = FALSE,
         annotation_col = ann,
         color = cols, border_color = NA,
         annotation_colors = ann_cols)
dev.off()

# figure_2B
png('manuscript/figures/figure2B.png', width = 3, height = 3, units = 'in', res = 500, pointsize = 6)
plotDispEsts(dds, ylab = 'Dispersion', xlab = 'Mean of normalized counts', cex.axis = .7)
dev.off()

# figure_2C
png('manuscript/figures/figure2C.png', width = 3, height = 3, units = 'in', res = 500, pointsize = 6)
sample_dist <- dist(t(assay(dds)))
mds <- data.frame(cmdscale(as.matrix(sample_dist)))
plot(mds$X1, mds$X2, pch = 19, col = sample_table$stage,
     xlab = 'M1', ylab = 'M2')
legend('bottomright', legend = unique(sample_table$stage), pch = 19, col = unique(sample_table$stage), cex = .8)
dev.off()

# figure_2D
png('manuscript/figures/figure2D.png', width = 3, height = 3, units = 'in', res = 500, pointsize = 6)
vdf <- res_groups %>%
  filter(stage == 'stage_Differentiation_vs_Proliferation') %>%
  mutate(data_set = 'Al Adhami et al') %>%
  bind_rows(res2) %>%
  filter(row %in% unlist(go[['GO:0006914']])) %>%
  dcast(row ~ data_set, value.var = 'baseMean')

plot(vdf$`Al Adhami et al`, vdf$`Duteil et al`,
     pch = 19,
     cex.axis = .7,
     xlab = 'Al Adhami et al', ylab = 'Duteil et al')
dev.off()

#figure_3
## figure: autophagy clusters
e <- assay(dds)
e <- (e - rowMeans(e))/rowSds(e)

df <- melt(e) %>%
  setNames(c('gene_id', 'sample', 'count')) %>%
  left_join(sample_table) %>%
  group_by(stage, gene_id) %>%
  summarise(count = mean(count)) %>%
  right_join(cm_sim) %>%
  group_by(cluster, stage) %>%
  mutate(ymin = min(count), ymax = max(count))

ggplot(df, aes(x = stage, y = count, group = gene_id)) +
  facet_wrap(~cluster, ncol = 4) +
  geom_ribbon(data = subset(df, cluster %in% paste(1:4)),
              aes(ymin = ymin, ymax = ymax),
              fill = 'gray', alpha = .5) +
  geom_smooth(data = subset(df, cluster %in% paste(1:4)),
              aes(group = cluster)) +
  geom_line(data = subset(df, cluster %in% paste('Similar', 1:4)), 
            aes(color = gene_id)) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 1),
        strip.background = element_rect(fill = 'gray', colour = 'white'),
        legend.position = 'none') +
  labs(x = '', y = 'Standarized Expression')
ggsave('manuscript/figures/figure3.png', width = 6, height = 4, units = 'in', dpi = 500)

# figure_4
res_go %>%
  filter(category %in% unlist(mget('GO:0006914', GOBPOFFSPRING))) %>%
  dplyr::select(stage, category, starts_with('Prop'), FDR) %>%
  gather(prop, val, PropDown, PropUp) %>%
  mutate(val = ifelse(prop == 'PropDown', -val, val),
         Signifaicant = ifelse(FDR < .1, 'Yes', 'NO'),
         stage = factor(stage, levels = unique(stage))) %>%
  ggplot(aes(x = category, y = val, fill = Signifaicant)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  facet_wrap(~stage, nrow = 3, strip.position = 'left') +
  theme_minimal() +
  theme(legend.position = 'top',
        legend.justification = 'right', 
        legend.text = element_text(size = 6),
        legend.key.size = unit(.3, units = 'cm'),
        legend.title = element_text(size = 6),
        axis.text.x = element_text(angle = 60, size = 7, hjust = 1),
        strip.background = element_rect(fill = 'lightgray', color = NA)) +
  labs(x = '', y = 'Proportion of DE Genes')
ggsave('manuscript/figures/figure4.png', width = 6.5, height = 5, units = 'in', dpi = 500)

# figure 5  
gene_data <- fpm(dds) %>%
  melt %>%
  setNames(c('gene_id', 'sample', 'fpm')) %>%
  left_join(sample_table) %>%
  group_by(gene_id, stage) %>%
  summarise(ave = mean(fpm)) %>%
  acast(gene_id ~ stage) %>%
  scale

select <- AnnotationDbi::select
pathview(gene_data[,c(3,4)], pathway.id = '04140', species = 'Mus musculus',
         multi.state = TRUE, gene.idtype = 'SYMBOL', out.suffix = 'autophagy23')

# figure S3
png('manuscript/figures/figureS3.png', width = 6.5, height = 6.5, units = 'in', res = 500)

pathways <- c('04150', '04630', '04910', '04920')
names(pathways) <- LETTERS[1:4]
par(mfrow=c(2,2), mar = c(4, 4, 2, 2))
for(i in seq_along(pathways)) {
  p <- get_pathway(res_groups, 'stageDifferentiation', pathways[i], kegg)
  barcodeplot(p$pvalues,
              index = p$index,
              gene.weights = p$weights,
              quantiles = c(.5, .1),
              labels = c("Low", "High"),
              weights.label = 't-stat',
              xlab = 'p-values')
  mtext(names(pathways)[i], line = 0, adj = -.2)
}
dev.off()

#figure S4
plot.new()
png('manuscript/figures/figureS4.png', width = 6.5, height = 6, units = 'in', res = 500)
plotDEXSeq(res_exon,
           'chr2_Acbd5+',
           fitExpToVar = 'stage',
           expression = FALSE,
           splicing = TRUE,
           displayTranscripts = TRUE,
           main = '',
           cex = .7, cex.lab = .7, cex.axis = .7)
mtext('Proliferation', line = 26, adj = 1, cex = .7, col = 'green')
mtext('Quiescence', line = 25, adj = 1, cex = .7, col = 'purple')
mtext('Clonal Expansion', line = 24, adj = 1, cex = .7, col = 'red')
mtext('Differentiation', line = 23, adj = 1, cex = .7,col = 'blue')
dev.off()
