#!/bin/bash

HANDLE="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"

# make directories for genome and data
mkdir mm10_iGenome
mkdir data

# get mm10 genome.fa and genes.gtf files
cd mm10_iGenome
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -xvzf Mus_musculus/UCSC/mm10/WholeGenomeFasta/genome.fa .
tar -xvzf Mus_musculus/UCSC/mm10/Annotations/Genes/genes.gtf .
cd ..

# get fastq.gz files for Al Adhami et al
cd data
for i in SRR969527 SRR969528 SRR969529 SRR969530 SRR969531 SRR969532 SRR969533 SRR969534
do
  wget -b $HANDLE/SRR969/$i/$i.fastq.gz &
done

# get fastq.gz files for Duteil et al
for i in SRR988305 SRR988306 SRR988307 SRR988308 SRR988309 SRR988310
do
  wget -b $HANDLE/SRR988/$i/$i_1.fastq.gz &
  wget -b $HANDLE/SRR988/$i/$i_2.fastq.gz &
done
cd ..