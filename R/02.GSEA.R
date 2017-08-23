# load data
load('processed_data.RData')

## figure: autophgy gene sets
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
        axis.text.x = element_text(angle = 90),
        strip.background = element_rect(fill = 'lightgray', color = NA)) +
  labs(x = '', y = 'Proportion of DE Genes (%)')

## KEGG map
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
         multi.state = TRUE, gene.idtype = 'SYMBOL',
         kegg.dir = 'output/kegg', out.suffix = 'autophagy23')

