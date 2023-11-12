# module load R/3.6.1
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")

# Rscript extract_spliceogen_and_spliceai_variants.R infile outfile
# Extract splicing variants, selecting by spliceogen score, spliceai score, spliceaiIndels score, and filter by gnomad.
# Also extract multi-nucleotide variants called by spliceogen score, for which we know we don't have a pre-computed spliceai score, filter by gnomad,
# to be sent as input to spliceai-tensorflow to get a spliceai score.
#
# Select splice site loss variants as any of the following:
# (spliceogen.donLossP >= 0.8) and (withinSite != ".")
# (spliceogen.accLossP >= 0.8) and (withinSite != ".") (make spliceogen accepter cut-off 0.8 instead of 0.85 so as to make things simpler)
# spliceai_donor_loss_delta_score >= 0.5 (spliceai_score >= 0.8 is very high confidence, but 0.5 to 0.8 has true positive without too many false positives)
# spliceai_acceptor_loss_delta_score >= 0.5 (spliceai_score >= 0.8 is very high confidence, but 0.5 to 0.8 has true positive without too many false positives)
# spliceaiIndels_donor_loss_delta_score >= 0.5 (spliceai_score >= 0.8 is very high confidence, but 0.5 to 0.8 has true positive without too many false positives)
# spliceaiIndels_acceptor_loss_delta_score >= 0.5 (spliceai_score >= 0.8 is very high confidence, but 0.5 to 0.8 has true positive without too many false positives)
#
# Select splice site gain variants as any of the following:
# (spliceai_acceptor_gain_delta_score >= 0.5)
# (spliceai_donor_gain_delta_score >= 0.5)
# (spliceaiIndels_acceptor_gain_delta_score >= 0.5)
# (spliceaiIndels_donor_gain_delta_score >= 0.5)

# (Don't use spliceogen.donGainP or spliceogen.accGainP to choose splice site gain variants because we would only choose them if spliceai_score was above threshold too)
#
# Steve's comments: That spec sounds good to me. It definitely won't catch all splicing variants. 
# When I have assessed datasets of known splice altering variants, many have lower SpliceAI scores- there will be many in the range of 0.2-0.5, and lower. 
# Of course, lowering the threshold is a trade-off of allowing many more false positive calls in. I think a threshold of 0.5 is a good option- the previous 0.8 was definitely too stringent.
# Actually, no need to filter for within splice site. The Spliceogen algorithm does not give any loss scores for variants if they are not within a splice site.
#
# For Spliceogen, due to the way it scans for motifs, when there are within-splice-site gain events, bascially the loss score is kind of meaningless.
# It's not uncommon for a variant to both disrupt an existing splice site and create a new splices site several BPs from the original.
# For variants within splice sites, the idea is to take the maximum of the Spliceogen gain and loss scores.
# Or just use the withinSS file, because that is what it already does for those columns. Either way is fine.
# Steve would consider the main score to just be the max of either score. In that case, it would indicate something atypical is going on, and we would need to take a look at the variant sequence manually.

# Also:
# Select splice site loss variants as any of the following:
# if (len(ref)==1)&(len(alt)>2) then is_non_spliceai_MNV = True
# if (len(ref)>5)&(len(alt)==1) then is_non_spliceai_MNV = True
# if (len(ref)>1)&(len(alt)>1) then is_non_spliceai_MNV = True
# is_non_spliceai_MNV and (spliceogen.donLossP >= 0.8) and (withinSite != ".")
# is_non_spliceai_MNV and (spliceogen.accLossP >= 0.8) and (withinSite != ".") (make spliceogen accepter cut-off 0.8 instead of 0.85 so as to make things simpler)
# is_non_spliceai_MNV and (spliceogen.donGainP >= 0.8)
# is_non_spliceai_MNV and (spliceogen.accGainP >= 0.8)

# Actually, when spliceogen calls a splice site loss variant due to donLossP or accLossP scores, then we consider it to be a true positive.
# However, when spliceogen calls a splice site gain variant due to donGainP or accGainP scores, 
# then we need to get spliceai score, from spliceai tensorflow if necessary, to confirm that it is a true positive.

options(width=200)
options(stringsAsFactors = FALSE)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( 'infile.txt', 'output_splicing_variants.txt', '0.05' ) # for testing
infile = as.character(args[1])
outfile = as.character(args[2])
gnomad_filter = as.numeric(args[3])

file_info = file.info(infile)
is_file_empty = file_info$size
if (is_file_empty == 0) {

  specific_input_is_empty = paste('Input file is empty: ',infile,sep='')
  input_is_empty = "CHROM\tPOS\tEND\tREF\tALT\tthere_are_no_variants"

  print( specific_input_is_empty )

  writeLines( input_is_empty, outfile )

  quit()
}

data = read.table( infile, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )

if (nrow(data) == 0) {

  specific_input_is_empty = paste('Input file has no data: ',infile,sep='')
  input_is_empty = "CHROM\tPOS\tEND\tREF\tALT\tthere_are_no_variants"

  print( specific_input_is_empty )

  writeLines( input_is_empty, outfile )

  quit()
}

#data[is.na(data)] = ''
#data[(data=='')] = '.'

# Rename any unreadable sample_format fields to readable names
names(data)[names(data) == 'FORMAT_AD'] = 'Allelic_depths_ref_alt'
names(data)[names(data) == 'FORMAT_DP'] = 'Approximate_read_depth'
names(data)[names(data) == 'FORMAT_GQ'] = 'Genotype_quality'
names(data)[names(data) == 'FORMAT_GT'] = 'Genotype'
names(data)[names(data) == 'FORMAT_MIN_DP'] = 'Min_depth_within_GVCF_block'
names(data)[names(data) == 'FORMAT_PGT'] = 'Physical_phasing_haplotype'
names(data)[names(data) == 'FORMAT_PID'] = 'Physical_phasing_ID'
names(data)[names(data) == 'FORMAT_PL'] = 'Phred_scaled_genotype_likelihoods'
names(data)[names(data) == 'FORMAT_RGQ'] = 'Reference_genotype_phred_confidence'
names(data)[names(data) == 'FORMAT_SB'] = 'Per_sample_component_statistics'

rename_colname_ensuring_it_is_unique = function( in_colnames, in_col ) {

	# G_AF.x => G_AF
	if (nchar(in_col) > 2) {
		if (str_sub(in_col,-2,-1) == '.x') {
			new_col = substr( in_col, 1, nchar(in_col)-2 )
			if (new_col %in% in_colnames) {
				dont_rename = 1
			} else {
				in_colnames[in_colnames == in_col] = new_col
				in_col = new_col
			}
		}
	}

	# INFO_G_AF => G_AF
	if (nchar(in_col) > 5) {
		if (substr(in_col, 1, 5) == 'INFO_') {
			new_col = substr( in_col, 6, nchar(in_col) )
			if (new_col %in% in_colnames) {
				dont_rename = 1
			} else {
				in_colnames[in_colnames == in_col] = new_col
				in_col = new_col
			}
		}
	}

	# G_AF.y => G_AF_repeated
	if (nchar(in_col) > 2) {
		if (str_sub(in_col,-2,-1) == '.y') {
			new_col = substr( in_col, 1, nchar(in_col)-2 )
			new_col = paste( new_col, '_repeated', sep='' )
			if (new_col %in% in_colnames) {
				dont_rename = 1
			} else {
				in_colnames[in_colnames == in_col] = new_col
				in_col = new_col
			}
		}
	}

	in_colnames
}

colnames_of_spliceogen_score_columns = c("spliceogen_donLossP", "spliceogen_accLossP", "spliceogen_donGainP", "spliceogen_accGainP")

colnames_of_spliceai_score_columns = c("spliceai_acceptor_gain_delta_score", "spliceai_donor_gain_delta_score", "spliceai_acceptor_loss_delta_score", "spliceai_donor_loss_delta_score", "spliceaiIndels_acceptor_gain_delta_score", "spliceaiIndels_donor_gain_delta_score", "spliceaiIndels_acceptor_loss_delta_score", "spliceaiIndels_donor_loss_delta_score", "sai_AG_DS", "sai_AL_DS", "sai_DG_DS", "sai_DL_DS", "saiI_AG_DS", "saiI_AL_DS", "saiI_DG_DS", "saiI_DL_DS", "saiM_AG_DS", "saiM_AL_DS", "saiM_DG_DS", "saiM_DL_DS", "saiIM_AG_DS", "saiIM_AL_DS", "saiIM_DG_DS", "saiIM_DL_DS")

colnames_of_mmsplice_score_columns = c("mmsplice_delta_logit_psi")

colnames_of_score_columns = c(colnames_of_spliceogen_score_columns, colnames_of_spliceai_score_columns, colnames_of_mmsplice_score_columns)

colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'gnomad211_VC_exome_AF' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'gnomad211_exome_AF' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'gnomad31_VC_genome_AF' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'gnomadMNVdist_AF_mnv_adj' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'spliceogen_withinSite' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'spliceogen_accGainP' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'spliceogen_donGainP' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'spliceogen_accLossP' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'spliceogen_donLossP' )
colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), 'spliceai_site' )

for (this_colname in colnames_of_score_columns) {
  colnames(data) = rename_colname_ensuring_it_is_unique( colnames(data), this_colname )
}

for (this_colname in colnames_of_score_columns) {
  if (this_colname %in% colnames(data)) {
    data[[this_colname]] = ifelse( is.na(data[[this_colname]]), "-1", data[[this_colname]] )
    data[[this_colname]] = ifelse( (data[[this_colname]]==""), "-1", data[[this_colname]] )
    data[[this_colname]] = ifelse( (data[[this_colname]]=="."), "-1", data[[this_colname]] )
    data[[this_colname]] = as.numeric(as.character(data[[this_colname]]))
  }
}

data$spliceogen_withinSite[is.na(data$spliceogen_withinSite)] = ""

# Filter by whether this variant is a multi-nucleotide variant for which there is not a precomputed spliceai score.

# Set colnames to REF and ALT
names(data)[names(data) == 'ref'] = 'REF'
names(data)[names(data) == 'Ref'] = 'REF'
names(data)[names(data) == 'alt'] = 'ALT'
names(data)[names(data) == 'Alt'] = 'ALT'

data$is_non_spliceai_MNV = 0
data$is_non_spliceai_MNV = ifelse( ((nchar(data$REF)==1) & (nchar(data$ALT)>2)), 1, data$is_non_spliceai_MNV )
data$is_non_spliceai_MNV = ifelse( ((nchar(data$REF)>5) & (nchar(data$ALT)==1)), 1, data$is_non_spliceai_MNV )
data$is_non_spliceai_MNV = ifelse( ((nchar(data$REF)>1) & (nchar(data$ALT)>1)), 1, data$is_non_spliceai_MNV )

# Filter by splicing scores

# data2 = data[ ( ((data$spliceogen_donLossP>=0.8)&(data$spliceogen_withinSite!=".")&(data$spliceogen_withinSite!="")) | ((data$spliceogen_accLossP>=0.8)&(data$spliceogen_withinSite!=".")&(data$spliceogen_withinSite!="")) | (data$spliceai_donor_loss_delta_score>=0.5) | (data$spliceai_acceptor_loss_delta_score>=0.5) | (data$spliceaiIndels_donor_loss_delta_score>=0.5) | (data$spliceaiIndels_acceptor_loss_delta_score>=0.5) | (data$spliceai_acceptor_gain_delta_score>=0.5) | (data$spliceai_donor_gain_delta_score>=0.5) | (data$spliceaiIndels_acceptor_gain_delta_score>=0.5) | (data$spliceaiIndels_donor_gain_delta_score>=0.5) ), ]

for (this_colname in colnames_of_score_columns) {
  if (this_colname %in% colnames(data)) {
    data[[this_colname]] = ifelse( is.na(data[[this_colname]]), "-1", data[[this_colname]] )
    data[[this_colname]] = ifelse( (data[[this_colname]]==""), "-1", data[[this_colname]] )
    data[[this_colname]] = ifelse( (data[[this_colname]]=="."), "-1", data[[this_colname]] )
    data[[this_colname]] = as.numeric(as.character(data[[this_colname]]))
  }
}

data$confident_splice_variant_scores = 0

for (this_colname in colnames_of_spliceogen_score_columns) {
  if (this_colname %in% colnames(data)) {
    data$confident_splice_variant_scores = ifelse( (data[[this_colname]] >= 0.8), 1, data$confident_splice_variant_scores )
  }
}

for (this_colname in colnames_of_spliceai_score_columns) {
  if (this_colname %in% colnames(data)) {
    data$confident_splice_variant_scores = ifelse( (data[[this_colname]] >= 0.5), 1, data$confident_splice_variant_scores )
  }
}

# Filter by whether spliceogen thinks that this multi-nucleotide variant is a splicing variant.
# data2 = data[ ((data$is_non_spliceai_MNV==1) & ( ((data$spliceogen_donLossP>=0.8)&(data$spliceogen_withinSite!=".")&(data$spliceogen_withinSite!="")) | ((data$spliceogen_accLossP>=0.8)&(data$spliceogen_withinSite!=".")&(data$spliceogen_withinSite!="")) | (data$spliceogen_donGainP>=0.8) | (data$spliceogen_accGainP>=0.8) )), ]

data$needs_spliceai_tensorflow_scores = 0

for (this_colname in colnames_of_spliceogen_score_columns) {
  if (this_colname %in% colnames(data)) {
    data$needs_spliceai_tensorflow_scores = ifelse( ((data$is_non_spliceai_MNV==1) & (data[[this_colname]] >= 0.8)), 1, data$needs_spliceai_tensorflow_scores )
  }
}

for (this_colname in colnames_of_mmsplice_score_columns) {
  if (this_colname %in% colnames(data)) {
    data$needs_spliceai_tensorflow_scores = ifelse( ( (data$is_non_spliceai_MNV==1) & ((data[[this_colname]] >= 2) | (data[[this_colname]] <= -2)) ), 1, data$needs_spliceai_tensorflow_scores )
  }
}

data2 = data[ ((data$confident_splice_variant_scores==1) | (data$needs_spliceai_tensorflow_scores==1)), ]

# Filter by gnomad frequencies

colnames_of_gnomad_frequency_columns = c("gnomad211_VC_exome_AF", "gnomad211_exome_AF", "gnomad31_VC_genome_AF", "gnomadMNVdist_AF_mnv_adj")

data2$filter_by_frequency = 0
data2$filter_by_frequency = as.numeric(data2$filter_by_frequency)
if ('gnomad211_VC_exome_AF' %in% colnames(data2)) {
	data2$gnomad211_VC_exome_AF[is.na(data2$gnomad211_VC_exome_AF)] = "-1"
	data2$gnomad211_VC_exome_AF = ifelse( (data2$gnomad211_VC_exome_AF==""), "-1", data2$gnomad211_VC_exome_AF )
	data2$gnomad211_VC_exome_AF = ifelse( (data2$gnomad211_VC_exome_AF=="."), "-1", data2$gnomad211_VC_exome_AF )
	data2$gnomad211_VC_exome_AF = as.numeric(as.character(data2$gnomad211_VC_exome_AF))
	data2$filter_by_frequency = apply( data2[,c("gnomad211_VC_exome_AF","filter_by_frequency")], 1, max )
}
if ('gnomad211_exome_AF' %in% colnames(data2)) {
	data2$gnomad211_exome_AF[is.na(data2$gnomad211_exome_AF)] = "-1"
	data2$gnomad211_exome_AF = ifelse( (data2$gnomad211_exome_AF==""), "-1", data2$gnomad211_exome_AF )
	data2$gnomad211_exome_AF = ifelse( (data2$gnomad211_exome_AF=="."), "-1", data2$gnomad211_exome_AF )
	data2$gnomad211_exome_AF = as.numeric(as.character(data2$gnomad211_exome_AF))
	data2$filter_by_frequency = apply( data2[,c("gnomad211_exome_AF","filter_by_frequency")], 1, max )
}
if ('gnomad31_VC_genome_AF' %in% colnames(data2)) {
	data2$gnomad31_VC_genome_AF[is.na(data2$gnomad31_VC_genome_AF)] = "-1"
	data2$gnomad31_VC_genome_AF = ifelse( (data2$gnomad31_VC_genome_AF==""), "-1", data2$gnomad31_VC_genome_AF )
	data2$gnomad31_VC_genome_AF = ifelse( (data2$gnomad31_VC_genome_AF=="."), "-1", data2$gnomad31_VC_genome_AF )
	data2$gnomad31_VC_genome_AF = as.numeric(as.character(data2$gnomad31_VC_genome_AF))
	data2$filter_by_frequency = apply( data2[,c("gnomad31_VC_genome_AF","filter_by_frequency")], 1, max )
}
if ('gnomadMNVdist_AF_mnv_adj' %in% colnames(data2)) {
        data2$gnomadMNVdist_AF_mnv_adj[is.na(data2$gnomadMNVdist_AF_mnv_adj)] = "-1"
        data2$gnomadMNVdist_AF_mnv_adj = ifelse( (data2$gnomadMNVdist_AF_mnv_adj==""), "-1", data2$gnomadMNVdist_AF_mnv_adj )
        data2$gnomadMNVdist_AF_mnv_adj = ifelse( (data2$gnomadMNVdist_AF_mnv_adj=="."), "-1", data2$gnomadMNVdist_AF_mnv_adj )
        data2$gnomadMNVdist_AF_mnv_adj = as.numeric(as.character(data2$gnomadMNVdist_AF_mnv_adj))
        data2$filter_by_frequency = apply( data2[,c("gnomadMNVdist_AF_mnv_adj","filter_by_frequency")], 1, max )
}
data3 = data2[ (data2$filter_by_frequency <= gnomad_filter), ]

# Sort from highest splice variant score to lower, then by lowest gnomad frequency to highest, which was filtered to be <= 5%.

#data3$sort_by_splice_score = apply( data3[,c("spliceogen_donLossP", "spliceogen_accLossP", "spliceogen_donGainP", "spliceogen_accGainP", "spliceai_acceptor_gain_delta_score", "spliceai_donor_gain_delta_score", "spliceai_acceptor_loss_delta_score", "spliceai_donor_loss_delta_score", "spliceaiIndels_acceptor_gain_delta_score", "spliceaiIndels_donor_gain_delta_score", "spliceaiIndels_acceptor_loss_delta_score", "spliceaiIndels_donor_loss_delta_score")], 1, max )

list_of_splicing_colnames = c()
for (this_colname in colnames_of_spliceogen_score_columns) {
  if (this_colname %in% colnames(data3)) {
    list_of_splicing_colnames = c(list_of_splicing_colnames, this_colname)
  }
}
for (this_colname in colnames_of_spliceai_score_columns) {
  if (this_colname %in% colnames(data3)) {
    list_of_splicing_colnames = c(list_of_splicing_colnames, this_colname)
  }
}

data3$sort_by_splice_score = apply( data3[, list_of_splicing_colnames ], 1, max )

data4 = data3[order( data3$filter_by_frequency, -data3$sort_by_splice_score ),]

data4$filter_by_frequency = NULL
data4$sort_by_splice_score = NULL

for (this_colname in colnames_of_gnomad_frequency_columns) {
  if (this_colname %in% colnames(data)) {
    data4[[this_colname]] = as.character(data4[[this_colname]])
    data4[[this_colname]] = ifelse( (data4[[this_colname]]=="-1"), "", data4[[this_colname]] )
  }
}

for (this_colname in colnames_of_score_columns) {
  if (this_colname %in% colnames(data)) {
    data4[[this_colname]] = as.character(data4[[this_colname]])
    data4[[this_colname]] = ifelse( (data4[[this_colname]]=="-1"), "", data4[[this_colname]] )
  }
}

data4[is.na(data4)] = ''
data4 = unique(data4)

write.table( data4, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


