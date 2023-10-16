#!/bin/bash

set -e

echo "Please enter your SRA:"
read SRR_id 

echo "Fetching fastq files"
fasterq-dump $SRR_id
echo "Input data collected successfully"

echo "Star is used to align with the reference genome"
STAR  --runThreadN 10 --genomeDir /biodb/human/gencode/v35/star_index --readFilesIn ${SRR_id}.fastq --outSAMtype BAM Unsorted --genomeLoad NoSharedMemory --limitOutSJcollapsed 5000000 --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts 

echo " Mudskipper is used to convert a bulk RNA-Seq genomic BAM to a transcriptomic BAM for quantification"
/data008/users/ttemker/Differences_in_quantification/data/mudskipper/target/release/mudskipper bulk --gtf /biodb/human/gencode/v35/gene_annotations.gtf --alignment Aligned.out.bam --out Aligned.mudskipperout.bam

echo " Starting to quantify mudskipper output with salmon"
salmon quant --libType A --alignments Aligned.mudskipperout.bam --targets /biodb/human/gencode/v35/transcripts.fa.gz --output salmon_mudskipper --threads 10 --gencode

echo "Starting to quantify STAR output with salmon"
salmon quant --libType A --alignments Aligned.toTranscriptome.out.bam --targets /biodb/human/gencode/v35/transcripts.fa.gz --output salmon_star --threads 10 --gencode 

echo "Starting salmon quant"
salmon quant --libType A --threads 10 --index /biodb/human/gencode/v35/salmon_index --output salmon_quant -r ${SRR_id}.fastq 

echo "done"

echo " Calculating correlation testing"
Rscript correlation.R
echo " Done "