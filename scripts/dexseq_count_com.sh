#!/bin/bash

SCRIPT="scripts/dexseq_count.py"
GTF="mm10_iGenome/genes.gtf"
GFF="mm10_iGenome/dexseq_genes.gff"
ALIGNMENT="output/hisat2_output"
OUTPUT="output/htseq_exon"

# flatten gtf file to use with dexseq-coun
python dexseq_prepare_annotation.py $GTF $GFF

# make directory of the count output
mkdir output/htseq_exon

# run dexseq for Al Adhami et al aligned reads
for i in SRR969527 SRR969528 SRR969529 SRR969530 SRR969531 SRR969532 SRR969533 SRR969534
do
	python $SCRIPT \
	-s no \
	$GFF \
	$ALIGNMENT/$i.sam \
	$OUTPUT/$i.txt
done
