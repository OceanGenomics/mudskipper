library(polyester)
library(Biostrings)

# read in the transcript fasta file
fasta = readDNAStringSet("/biodb/human/gencode/v35/transcripts.fa.gz")

# read in the quant files
quant72 = read.delim('/data102/users/gus/Downloads/quants/advanced_melanoma/sample_72/salmon/quant.sf')
quant73 = read.delim('/data102/users/gus/Downloads/quants/advanced_melanoma/sample_73/salmon/quant.sf')
# calculate I_bar
I_bar72 = sum(quant72$TPM/1000000*quant72$EffectiveLength)
I_bar73 = sum(quant73$TPM/1000000*quant73$EffectiveLength)
# calculate FPKM
FPKM72 = 1000/I_bar72*quant72$TPM
FPKM73 = 1000/I_bar73*quant73$TPM

# create FPKM matrix
fpkmDF = data.frame(FPKM72, FPKM73)
fpkmMat = data.matrix(fpkmDF)
rownames(fpkmMat) = quant72$Name

# change names for transcript fasta
new_names = unlist(lapply(names(fasta), function(i) {strsplit(i, "[|]")[[1]][1]}))
names(fasta) = new_names
# filter
# get the transcript names that exist in transcript.fa but not in quant.sf
delete_names = setdiff(names(fasta), quant72$Name)
new_transcript = fasta[!names(fasta) %in% delete_names]
writeXStringSet(new_transcript, "transcripts2.fa")

# run simulate_experiment_empirical
simulate_experiment_empirical(fpkmMat=fpkmMat, fasta="transcripts2.fa", grouplabels=c("sample72", "sample73"), outdir='simulated_reads_72_73')
