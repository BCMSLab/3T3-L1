#!/bin/bash

ALIGNMENT="output/hisat2_output"
GTF="mm10_iGenome/genes.gtf"
OUTPUT="output/htseq_gene"

# make directory of the count output
mkdir output/htseq_gene

# run htseq count for Al Adhami et al aligned reads
for i in SRR969527 SRR969528 SRR969529 SRR969530 SRR969531 SRR969532 SRR969533 SRR969534
do
	htseq-count -s no \
	$ALIGNMENT/$i.sam \
	$GTF \
	> $OUTPUT/$i.txt
done

# run htseq count for Duteil et al aligned reads
for i in SRR988305 SRR988306 SRR988307 SRR988308 SRR988309 SRR988310
do
	htseq-count -s no \
	htseq-count -s no \
        $ALIGNMENT/$i.sam \
        $GTF \
        > $OUTPUT/$i.txt
done

