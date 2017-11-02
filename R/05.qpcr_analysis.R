library(readr)
library(dplyr)
library(pcr)
library(ggplot2)

ct <- read_csv('data/qpcr_ct.csv')
names(ct)[1:2] <- c('S_18', 'B_actine')

group <- rep(c('Proliferation', 'Quiescence', 'Clonal Expansion', 'Differentiation'), each = 3)

dat <- pcr_analyze(ct, group_var = group,
                   reference_gene = 'S_18',
                   reference_group = 'Proliferation')
dat %>%
  filter(gene != 'B_actine') %>%
  mutate(group = factor(group, levels = c('Proliferation', 'Quiescence', 'Clonal Expansion', 'Differentiation')),
         gene = factor(gene, levels = c('Ubqln2', 'Pink1', 'Foxo1', 'Prkaa2'))) %>%
  ggplot(aes(x = group, y = norm_rel, group = gene)) +
  geom_line() +
  facet_wrap(~gene, nrow = 2) +
  lims(y = c(0, 4.2)) +
  geom_errorbar(aes(ymin = int_lower, ymax = int_upper),  width = .5) +
    theme_bw()+
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.title.y = element_text(size = 8),
          plot.margin = unit(c(0,0,0,.5), 'cm'),
          panel.grid = element_blank(),
          panel.border = element_rect(size = 1)) +
    labs(x = '', y = 'Relative mRNA expression')

ggsave('manuscript/figures/new_fig4.png', width = 6, height = 6, units = 'in', dpi = 500)

dat <- pcr_caliberate(ct[, 1:2], group_var = group, reference_group = 'stage_1')
dat %>%
  ggplot(aes(x = group, y = rel)) +
  geom_col() +
  facet_wrap(~gene) +
  geom_errorbar(aes(ymin = int_lower, ymax = int_upper))
