# read in the quant files
quant_a = read.delim('salmon_star/quant.sf')
quant_m = read.delim('salmon_mudskipper/quant.sf')
# filter the transcripts
# keep transcripts exist in both files
delete_names = setdiff(quant_a$Name, quant_m$Name)
new_quant_a = quant_a[!quant_a$Name %in% delete_names,]
corr = cor(new_quant_a$TPM, quant_m$TPM, method = "spearman")
cat("Correlation value of star vs mudskipper:", corr, "\n")
cat("\014") 

# read in the quant files
quant_x = read.delim('salmon_mudskipper/quant.sf')
quant_y = read.delim('salmon_star/quant.sf')
# filter the transcripts
# keep transcripts exist in both files
delete_names_a = setdiff(quant_x$Name, quant_y$Name)
new_quant_b = quant_x[!quant_x$Name %in% delete_names,]
corr_a = cor(new_quant_b$TPM, quant_y$TPM, method = "spearman")
cat("Correlation value of mudskipper vs star:", corr_a, "\n")
cat("\014") 

# read in the quant files
quant_c = read.delim('salmon_mudskipper/quant.sf')
quant_d = read.delim('salmon_quant/quant.sf')

# filter the transcripts
# keep transcripts exist in both files
delete_names_b = setdiff(quant_c$Name, quant_d$Name)
new_quant_c = quant_c[!quant_c$Name %in% delete_names_b,]

corr_b = cor(new_quant_c$TPM, quant_d$TPM, method = "spearman")
cat("Correlation value of mudskipper vs quant:", corr_b, "\n")
