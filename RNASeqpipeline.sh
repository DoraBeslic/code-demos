#!/bin/bash 
# script to generate a count matrix from bulk RNA-Seq reads

# demo data (bulk RNA-Seq reads):
# https://drive.google.com/file/d/1DGHjbhcRy_zTm6H9C_AUpkzBML-JhtA3/view

# change working directory
mkdir RNASeq_pipeline
cd ~/Desktop/code_demos/RNASeq_pipeline


# STEP 1: Quality control and pre-processing------------------------------------

# run fastqc to check quality of fastq reads
fastqc data/demo.fastq -o data/

# run Trimmomatic to trim reads with poor quality
java -jar /Users/isidorabeslic/Desktop/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 data/demo.fastq data/demo_trimmed.fastq TRAILING:10 -phred33
echo "Trimmomatic finished running!"

# check quality after trimming
fastqc data/demo_trimmed.fastq -o data/


# STEP 2: Alignment with HISAT2-------------------------------------------------

# get the genome indices
mkdir HISAT2
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -zxvf grch38_genome.tar.gz -C HISAT2/
rm grch38_genome.tar.gz

# run alignment
hisat2 -q --rna-strandness R -x HISAT2/grch38/genome -U data/demo_trimmed.fastq | samtools sort -o HISAT2/demo_trimmed.bam
echo "HISAT2 finished running!"


# STEP 3: Quantification with featureCounts-------------------------------------

# get gtf
wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
gzip -d Homo_sapiens.GRCh38.106.gtf.gz

# run quantification
mkdir quants
featureCounts -a Homo_sapiens.GRCh38.106.gtf -o quants/demo_featurecounts.txt HISAT2/demo_trimmed.bam
echo "FeatureCounts finished running!"
