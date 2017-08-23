# laod data
load('processed_data.RData')

## figure: numbers of differentially expressed genes
ind <- gene_annotation %>%
  filter(GO == 'GO:0006914') %>%
  dplyr::select(SYMBOL) %>%
  unlist %>%
  unique

cols <- c('#FEE6CE', '#BDD7E7', '#65BC6E')
res_groups %>%
  filter(row %in% ind) %>%
  mutate(stage = str_sub(stage, start = 6),
         padj = ifelse(padj < .1, 1, 0)) %>%
  acast(row ~ stage, value.var = 'padj') %>%
  modified_vennDiagram(circle.col = c('red', 'blue', 'green'),
                       mar = rep(0, 4),
                       cex = rep(.7, 4))

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
                  `Time Point` = sample_table$lab,
                  row.names = sample_table$sample)
ann_cols <- colorRampPalette(brewer.pal(n = 4, name = 'Blues'))(4) 
names(ann_cols) <- unique(ann$`Time Point`)
ann_cols <- list(Induction = c(Non = '#FEE6CE', IDX = '#E6550D'),
                 `Time Point` = c(`- 0 Days` = '#EFF3FF',
                                  `- 2 Days` = '#BDD7E7',
                                  `10 Hours` = '#6BAED6',
                                  `6 Days` = '#2171B5'))
pheatmap(mat,
         scale = 'row', legend = FALSE,
         cellheight = 15, cellwidth = 20, fontsize = 10, cluster_cols = FALSE,
         annotation_col = ann,
         color = cols, border_color = NA,
         annotation_colors = ann_cols)

## figure: EDA, model daignostics and validation
## A) Dispersion
par(mfrow=c(1,3))

plotDispEsts(dds)

## B) MDS
sample_dist <- dist(t(assay(dds)))
mds <- data.frame(cmdscale(as.matrix(sample_dist)))
plot(mds$X1, mds$X2, pch = 19, col = sample_table$stage,
     xlab = 'M1', ylab = 'M2')
legend('bottomright', legend = unique(sample_table$stage), pch = 19, col = unique(sample_table$stage))

## C) validataion
vdf <- res_groups %>%
  filter(stage == 'stageDifferentiation') %>%
  dplyr::select(-stage) %>%
  mutate(data_set = 'Adhami et al') %>%
  bind_rows(res2) %>%
  filter(row %in% ind) %>%
  dcast(row~data_set, value.var = 'log2FoldChange')
lims <- c(-3, 3)
plot(vdf$`Adhami et al`, vdf$`Duteil et al`, pch = 19,
     xlab = 'Adhami et al', ylab = 'Duteil et al',
     xlim = lims, ylim = lims)
abline(0, 1, col = 'gray')


## clustering
e <- assay(dds)[ind,]
e <- e[rowSds(e) > 100,]
e <- (e - rowMeans(e))/rowSds(e)

set.seed(123)
cm <- cmeans(e, 4, m = 1.5, dist = 'manhattan')
clusters <- data.frame(cluster = cm$cluster) %>%
  rownames_to_column('gene_id')

## figure: autophagy clusters
a <- as.data.frame(e) %>%
  rownames_to_column('gene_id') %>%
  gather(sample, exp, -gene_id) %>%
  left_join(clusters) %>%
  left_join(sample_table) %>%
  group_by(gene_id, stage) %>% mutate(exp = mean(exp)) %>%
  group_by(cluster, stage) %>% mutate(ymin = min(exp), ymax = max(exp)) %>%
  ggplot(aes(x = lab2, y = exp, ymin = ymin, ymax = ymax, group = cluster)) +
  geom_ribbon(fill = 'gray', alpha = .5) +
  geom_smooth() +
  facet_wrap(~cluster, ncol = 4) + 
  theme_minimal() +
  theme(strip.background = element_rect(fill = 'gray', colour = 'white'),
        axis.text.x = element_blank()) +
  labs(y = 'Standarized Expression', x = '')

## figure: similar genes
e <- assay(dds)
e <- e[rowSds(e) > 100,]
e <- (e - rowMeans(e))/rowSds(e)

mediods <- apply(cm$membership, 2, which.max)
mediods <- rownames(cm$membership)[mediods]

exc <- setdiff(ind, mediods)
e <- e[!rownames(e) %in% exc,]
ilst <- which(rownames(e) %in% mediods)
set.seed(123)
sim <- genefinder(e, ilist = ilst, numResults = 5)
sim <- melt(sim, value.name = 'index') %>%
  filter(L2 == 'indices') %>%
  mutate(gene_id = rownames(e)[index],
         cluster = L1) %>%
  dplyr::select(gene_id, cluster)

b <- as.data.frame(e) %>%
  rownames_to_column('gene_id') %>%
  gather(sample, exp, -gene_id) %>%
  right_join(sim) %>%
  left_join(sample_table) %>%
  group_by(gene_id, stage) %>%
  mutate(exp = mean(exp), cluster = as.character(cluster)) %>%
  ggplot(aes(x = lab2, y = exp, group = gene_id, color = gene_id)) +
  geom_line() +
  facet_wrap(~cluster, ncol = 4) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = 'gray', colour = 'white'),
        legend.position = 'none',
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  labs(y = 'Standarized Expression', x = '')
plot_grid(a, b, nrow = 2,
          labels = 'AUTO')
