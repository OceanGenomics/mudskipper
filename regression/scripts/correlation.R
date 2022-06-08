# read in the quant files
quant_a = read.delim('salmon_quant_a/quant.sf')
quant_m = read.delim('salmon_quant_m/quant.sf')

# filter the transcripts
# keep transcripts exist in both files
delete_names = setdiff(quant_a$Name, quant_m$Name)
new_quant_a = quant_a[!quant_a$Name %in% delete_names,]

corr = cor(new_quant_a$TPM, quant_m$TPM, method = "spearman")
if (corr < 0.90){
    quit(status=1)
}
