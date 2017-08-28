library(tidyverse)
library(reshape2)
list.files('output/htseq_gene/')
ll <- list()
ll$adhami <- data.frame(file = c('SRR969529.txt', 'SRR969530.txt', 'SRR969533.txt', 'SRR969534.txt'),
                     time = rep(c('proliferating', 'differentiatinh'), each = 2))
ll$duteil <- data.frame(file = c("SRR988305.txt", "SRR988306.txt", "SRR988307.txt", "SRR988308.txt", "SRR988309.txt", "SRR988310.txt"),
                     time = rep(c('proliferating', 'differentiatinh'), each = 3))
pdata <- bind_rows(ll, .id = 'dataset')

dat <- lapply(pdata$file, function(x) {
  read_tsv(paste('output/htseq_gene/', x, sep = ''), col_names = c('gene', 'count'))
})

names(dat) <- pdata$file
dat <- bind_rows(dat, .id = 'file')
head(dat)

gene_id <- c('Atg7', 'Atg5', 'Cebpa', 'Pparg', 'Rab7', 'Adrb2', 'Pink1', 'Park2', 'Fam134b', 'Optn', 'Maplc3a')

dat %>%
  filter(gene %in% gene_id) %>%
  left_join(pdata) %>%
  group_by(dataset, time, gene) %>%
  summarise(ave = log(mean(count)) + 1) %>%
  ggplot(aes(x = gene, y = ave, col = time)) +
  geom_point(position = 'dodge') +
  facet_wrap(~dataset, ncol = 1)

dat %>%
  filter(gene %in% gene_id) %>%
  left_join(pdata) %>%
  group_by(dataset, gene) %>%
  mutate(ave = scale(count)) %>%
  ggplot(aes(x = time, y = gene, fill = ave)) +
  geom_tile()+ facet_wrap(~dataset, ncol = 1)
