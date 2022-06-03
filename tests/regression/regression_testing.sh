#!/bin/bash

# Fail on error
set -e

# XXX that will need to be fixed. Hard coded paths
# add path for STAR and salmon
PATH=$PATH:/data102/users/gus/buildocean/local/bin
# add path for mudskipper
PATH=$PATH:/home/yfei/mudskipper/target/release

SCRIPT_DIR=$(dirname $(realpath -e ${BASH_SOURCE[0]:-$0}))

# XXX simulate_reads also has too many hard coded paths
# generate simulated reads
echo Simulate reads
Rscript ${SCRIPT_DIR}/scripts/simulate_reads.R


run STAR -> mudskipper -> Salmon pipeline
STAR
echo STAR
mkdir -p STAR_output
(cd STAR_output
 STAR --runThreadN 8 --genomeDir /biodb/human/gencode/v35/star_index \
     --readFilesIn ../simulated_reads_72_73/sample_01_1.fasta,../simulated_reads_72_73/sample_02_1.fasta ../simulated_reads_72_73/sample_01_2.fasta,../simulated_reads_72_73/sample_02_2.fasta
 )

# mudskipper
echo mudskipper
mudskipper bulk --gtf /biodb/human/gencode/v35/gene_annotations.gtf --alignment ./STAR_output/Aligned.out.sam --out transcriptomic.bam

# Salmon
echo salmon alignment
salmon quant \-t /biodb/human/gencode/v35/transcripts.fa.gz -l A -a transcriptomic.bam -o salmon_quant_a --gencode

# run Salmon with mapping mod
echo salmon mapping
salmon quant -i /biodb/human/gencode/v35/salmon_index -l IU --validateMappings -o salmon_quant_m \
-1 <(cat simulated_reads_72_73/sample_01_1.fasta simulated_reads_72_73/sample_02_1.fasta) \
-2 <(cat simulated_reads_72_73/sample_01_2.fasta simulated_reads_72_73/sample_02_2.fasta)

# calculate correlation
echo correlation
Rscript $SCRIPT_DIR/scripts/correlation.R
