# load data
load('processed_data.RData')

# numbers of differential exon usage
ind <- gene_annotation %>%
  filter(GO == 'GO:0006914') %>%
  dplyr::select(SYMBOL) %>%
  unlist %>%
  unique

df <- res_exon %>%
  as.data.frame %>%
  mutate(id = strsplit2(groupID, '_')[,2],
         id = str_sub(id, end = -2)) %>%
  filter(id %in% ind, padj < .1) %>%
  dplyr::select(groupID, id) %>%
  unique

res_exon %>%
  as.data.frame %>%
  filter(groupID %in% df$groupID) %>%
  dplyr::select(starts_with('log2fold')) %>%
  setNames(c('Quiescence', 'Clonal.expansion', 'Differentiation')) %>%
  mutate_all(function(x) ifelse(x > 1, 1, 0)) %>%
  vennDiagram(circle.col = c('red', 'blue', 'green'))


# plot of individual genes
plotDEXSeq(res_exon, 'chr1_Atg16l1+', fitExpToVar = 'stage',
           names = TRUE)

View(df)
df <- res_exon %>%
  as.data.frame %>%
  filter(groupID == 'chr1_Atg16l1+')
