#!/bin/bash

GENOME="mm10_iGenome/genom.fa"
INDEX='mm10_iGenome/hisat2_mm10_index'
DATA='data'
OUTPUT='output/hisat2_output'
P='4'

# build hista2 index
hisat2-build -p $P $GENOME $INDEX

# make directory of the alignment output
mkdir output/histat2_output

# run hisat2 for Al Adhami et al. dataset
for i in SRR969527 SRR969528 SRR969529 SRR969530 SRR969531 SRR969532 SRR969533 SRR969534
do 
	hisat2 -p $P \
	-q \
	-x $INDEX \
	-U $DATA/$i.fastq.gz \
	-S $OUTPUT/$i.sam
done

# run hisat2 for Duteil et al dataset
for i in SRR988305 SRR988306 SRR988307 SRR988308 SRR988309 SRR988310
do
	hisat2 -p $P \
	-q \
	-x $INDEX \
	-1 $DATA/$i_1.fastq.gz \
	-2 $DATA/$i_2.fastq.gz \
	-S $OUTPUT/$i.sam
done
