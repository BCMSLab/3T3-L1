library(DiagrammeR)
grViz("
      digraph main{
# nodes
## main
node [shape = box]
n1 [label = 'Raw Reads \n (ENA)']
n2 [label = 'Aligned reads \n (sam)']
n3 [label = 'Exon counts']
n4 [label = 'Gene counts']
n5 [label = 'Differential Expression']
n6 [label = 'Differential Exon \n Usage']
n7 [label = 'GSEA']
n8 [label = 'KEGG Pathways']
n9 [label = 'Autophagy \n gene sets']
n10 [label = 'Expression \n temporal \n profiles']
n11 [label = 'Genes with \n similar expression \n patterns']

## annotation
node [shape = box color = red]
a1 [label = 'USCS mm10 \n (genome.fa)']
a2 [label = 'mm10 index']
a3 [label = 'USCS mm10 \n (genes.gtf)']
a4 [label = 'flatten gtf \n (genes.gff)']

# connections
node [shape = none width = 0 height = 0 label = '']
p1
p2; p3; p4
p5; p6; p7

# ranks
{rank = same; n1; a1}
{rank = same; p1; a2}
{rank = same; n3; n4}
{rank = same; n5; n6}
{rank = same; n8; n9}
{rank = same; n10; n11}

{rank = same; p2; p3; p4}
{rank = same; p5; p6; p7}

{rankdir = 'RL'; a3; a4}
# edges
n1 -> p1 [arrowhead = none]
p1 -> n2 [label = 'HISAT2']
a1 -> a2 [label= 'hisat2-build']
a2 -> p1 [arrowhead = none minlen = 3]
n2 -> p3 [arrowhead = none]
p2 -> p3 -> p4 [arrowhead = none minlen = 2]
p2 -> n3 [label = 'dexseq-count']
n3 -> n6 [label = DEXSeq]
p4 -> n4 [label = 'HTSeq-count']
n4 -> n5 [label = DESeq2]
n5 -> n7 [label = limma]
p5 -> p6 -> p7 [arrowhead = none minlen = 2]
n7 -> p6 [arrowhead = none]
p5 -> n8
p7 -> n9
n9 -> n10 [label = 'C-means & genefilter']
n10 -> n11
      }
      ")
