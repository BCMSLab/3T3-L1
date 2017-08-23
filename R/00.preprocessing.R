################################################################################################################
# File Name: 00.preprocessing.R
# Description: This script does the preprocessing of the read counts of htseq and needed annotaions. Then,
#               performs the differential expression, exon usage and GESA.
# Created by: Mahmoud Ahmed <mahmoud.s.fahmy@students.kasralainy.edu.eg>
# Depends & input: output/htseq_gene/SRR9695*.txt, output/htseq_exon/SRR9695*.txt, mm10_iGenome/dexseq_genes.gtf
# Usage: Used with 00.functions.R & 00.figures.R to perform the analysis and make figures for T3-L1_manuscript
################################################################################################################

# load libraries
## libraries for analysis
library(DEXSeq)
library(DESeq2)
library(goseq)
library(limma)
library(genefilter)
library(Mfuzz)

## annotations libraries
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(GO.db)
library(pathview)

## libraries for data manipulation and plotting
library(tidyverse)
library(reshape2)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(DiagrammeR)
library(stringr)

## source functions
source('R/00.functions.R')

# Prepare gene annotation
gene_annotation <- AnnotationDbi::select(org.Mm.eg.db, keys(org.Mm.eg.db),
                                         c('SYMBOL', 'GENENAME', 'GO', 'PATH'),
                                         'ENTREZID')
go <- with(gene_annotation, split(SYMBOL, GO))
kegg <- with(gene_annotation, split(SYMBOL, PATH))

# differential expression
## Liklihood Ration Test (LRT)
sample_table <- data.frame(
    sample = strsplit2(list.files('output/htseq_gene', pattern = 'SRR9695*'), '\\.')[, 1],
    file_name = list.files('output/htseq_gene', pattern = 'SRR9695*'),
    stage = factor(rep(c('Proliferation', 'Quiescence', 'Clonal expansion', 'Differentiation'), each = 2),
                   levels = c('Proliferation', 'Quiescence', 'Clonal expansion', 'Differentiation')),
    induction = rep(c('Non', 'MDI'), each = 4)
    )

dds <- DESeqDataSetFromHTSeqCount(sample_table,
                                  directory = 'output/htseq_gene/',
                                  design = ~ stage)

dds <- dds[rowSums(counts(dds)) > 1,]
dds_groups <- DESeq(dds)
dds <- DESeq(dds,  test = 'LRT', reduced = ~1)
res_gene <- results(dds, tidy = TRUE)

## pairwise comparison
res_groups <- list()
for(i in 2:4) {
  res_groups[[i]] <- results(dds_groups,
                             name = resultsNames(dds_groups)[i],
                             tidy = TRUE)
  res_groups[[i]]$stage <- resultsNames(dds_groups)[i]
}
res_groups <- res_groups %>% bind_rows 


## validation data set (Duteil et al)
sample_table2 <-data.frame(sample = strsplit2(list.files('output/htseq_gene', pattern = 'SRR9883*'), '\\.')[, 1],
                           file_name = list.files('output/htseq_gene', pattern = 'SRR9883*'),
                           stage = rep(c('Day 0', 'Day 7'), each = 3))

dds2 <- DESeqDataSetFromHTSeqCount(sample_table2,
                                   directory = 'output/htseq_gene',
                                   design = ~stage)

dds2 <- dds2[rowSums(counts(dds2)) > 1,]
dds2 <- DESeq(dds2)
res2 <- results(dds2, tidy = TRUE)
res2$data_set <- 'Duteil et al'

# differential exon usage
## run dexseq
count_fls <- list.files('output/htseq_exon', pattern = 'SRR9695*', full.names = TRUE)

dxd <- DEXSeqDataSetFromHTSeq(count_fls,
                              sampleData = sample_table,
                              flattenedfile = 'mm10_iGenome/dexseq_genes.gff',
                              design = ~ sample + exon + stage:exon)

BPPARAM = MulticoreParam(workers = 4)
res_exon <- DEXSeq(dxd,
              fullModel = ~ sample + exon + stage:exon,
              reducedModel = ~ sample + exon,
              BPPARAM = BPPARAM,
              fitExpToVar = 'stage',
              quiet = TRUE)

ind <- res_exon %>%
  as.data.frame %>%
  mutate(SYMBOL = str_sub(strsplit2(groupID, '_')[,2], end = -2)) %>%
  filter(SYMBOL %in% unique(unlist(go[['GO:0006914']]))) %>%
  dplyr::select(groupID) %>%
  unlist %>%
  unique

## pairwise comparison
res_exon_groups <- list()
res_exon_groups$Quiescence <- dexseq_group(dxd, c(1:4, 9:12), 'Quiescence')
res_exon_groups$`Clonal expansion` <- dexseq_group(dxd, c(1:2, 5:6, 9:10, 13:14), 'Clonal expansion')
res_exon_groups$Differentiation <- dexseq_group(dxd, c(1:2, 7:8, 9:10, 15:16), 'Differentiation')

# GESEA
## gene ontology
model_matrix <- model.matrix(~sample_table$stage)
fpm <- fpm(dds)

go_mroast <- list()
for(i in 2:4){
  go_mroast[[i]] <- mroast(fpm, index = go, design = model_matrix, contrast = i)
  go_mroast[[i]]$stage <- levels(sample_table$stage)[i]
  go_mroast[[i]]$category <- rownames(go_mroast[[i]])
}

res_go <- go_mroast %>%
  bind_rows

## kegg analysis
kegg_mroast <- list()
for(i in 2:4) {
  kegg_mroast[[i]] <- mroast(fpm, index = kegg, design = model_matrix, contrast = i)
  kegg_mroast[[i]]$stage <- levels(sample_table$stage)[i]
  kegg_mroast[[i]]$category <- rownames(kegg_mroast[[i]])
}
res_kegg <- kegg_mroast %>%
  bind_rows

# clustering
## clustering autophagy genes
ind <- gene_annotation %>%
  filter(GO == 'GO:0006914') %>%
  dplyr::select(SYMBOL) %>%
  left_join(res_groups, by = c('SYMBOL'='row')) %>%
  filter(padj < .1, abs(log2FoldChange) > 1) %>%
  dplyr::select(SYMBOL) %>%
  unlist %>% unique

e <- assay(dds)[ind,]
e <- e[rowSds(e) > 100,]
e <- (e - rowMeans(e))/rowSds(e)

set.seed(123)
cm <- cmeans(e, 4, m = 1.5, dist = 'manhattan')
clusters <- data.frame(cluster = cm$cluster) %>%
  rownames_to_column('gene_id') %>%
  mutate(cluster = as.character(cluster))

## finding genes with similar expression patterns
e <- assay(dds)
e <- e[rowSds(e) > 100,]
e <- (e - rowMeans(e))/rowSds(e)

mediods <- apply(cm$membership, 2, which.max)
mediods <- rownames(cm$membership)[mediods]

exc <- setdiff(ind, mediods)
e <- e[!rownames(e) %in% exc,]
ilst <- match( mediods, rownames(e))

set.seed(123)
sim <- genefinder(e, ilist = ilst, numResults = 5)
names(sim) <- paste('Similar', c(1, 4, 2, 3))

sim <- melt(sim, value.name = 'index') %>%
  filter(L2 == 'indices') %>%
  mutate(gene_id = rownames(e)[index],
         cluster = L1) %>%
  dplyr::select(gene_id, cluster)

cm_sim <- bind_rows(clusters, sim)

# clean and save data
save.image('processed_data.RData')
