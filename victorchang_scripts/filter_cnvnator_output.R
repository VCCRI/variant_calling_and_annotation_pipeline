

options(width=200)
options(stringsAsFactors = FALSE)
library(tidyr)
library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( 'infile.tsv', 'outfile.tsv' ) # for testing

infile = as.character(args[1])
outfile = as.character(args[2])

data = read.table( infile, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )

colnames(data) = c("cnv_type", "cnv_coords", "cnv_size", "normalized_read_depth", "eval1_ttest", "eval2_gauss_prob", "eval3_ttest_midCNV", "eval4_gauss_prob_midCNV", "fraction_reads_mapped_with_q0_quality")

data$eval1_ttest = as.numeric(data$eval1_ttest)
data$eval2_gauss_prob = as.numeric(data$eval2_gauss_prob)
data$eval3_ttest_midCNV = as.numeric(data$eval3_ttest_midCNV)
data$eval4_gauss_prob_midCNV = as.numeric(data$eval4_gauss_prob_midCNV)

# It looks like this e-value means the same thing as in BLAST statistics: 
# the number of times we expect a hit of this significance would be observed by chance in a genome or database of this size. 
# For small values (e.g. below 0.05) the e-value and p-value converge on the same number, 
# but for p-values that approach 1.0, the e-values instead grow above 1. 
# An e-value >>1 means something similar to a p-value with leading 9's, i.e. almost certainly due to chance and not significant under the null model.

data2 = data[ ((data$eval1_ttest<=0.05)|(data$eval1_ttest>1)) & ((data$eval2_gauss_prob<=0.05)|(data$eval2_gauss_prob>1)) & ((data$eval3_ttest_midCNV<=0.05)|(data$eval3_ttest_midCNV>1)) & ((data$eval4_gauss_prob_midCNV<=0.05)|(data$eval4_gauss_prob_midCNV>1)),]

write.table( data2, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE )




