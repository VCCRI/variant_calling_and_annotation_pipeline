# module load R/3.6.1
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")

# Rscript extract_pathogenic_and_by_gene_name.R infile outfile genes_file

options(width=200)
options(stringsAsFactors = FALSE)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

args = commandArgs(trailingOnly=TRUE) # for production
#args = c('temp_60000_to_90000.tsv', 'temp_out_60000_to_90000.tsv', 'list_genes_CHD_all2021marV2_808_genes_withoutBars.txt') # for testing

infile = as.character(args[1])
outfile = as.character(args[2])
genes_file = as.character(args[3])

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
  write.table( data, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
  quit()
}

#data[is.na(data)] = ''
data[(data=='')] = '.'

# Remove the INFO_ prefix so that this code will work on data that doesn't have the INFO_ prefix.
# Make sure that we can tell the difference between VCF END and VCF INFO_END or any other INFO_ columns that might end up with the same name as some other column when we remove the INFO_
names(data)[names(data) == 'INFO_CHROM'] = 'VCF_INFO_CHROM'
names(data)[names(data) == 'INFO_POS'] = 'VCF_INFO_POS'
names(data)[names(data) == 'INFO_END'] = 'VCF_INFO_END'
names(data)[names(data) == 'INFO_ID'] = 'VCF_INFO_ID'
names(data)[names(data) == 'INFO_REF'] = 'VCF_INFO_REF'
names(data)[names(data) == 'INFO_ALT'] = 'VCF_INFO_ALT'
names(data)[names(data) == 'INFO_QUAL'] = 'VCF_INFO_QUAL'
names(data)[names(data) == 'INFO_FILTER'] = 'VCF_INFO_FILTER'
colnames(data) = gsub('^INFO_', '', colnames(data))

clinvar_date = ''
# look for eg. 'CLINVAR20200629_GENEINFO' in colnames. If multiples found then the last one found will be used.
# Rename the latest clinvar columns to remove date so that this and other programs know what the column names are, regardless of clinvar date.
for (this_colname in colnames(data)) {
  if ((str_sub(this_colname, 1, 7) == "CLINVAR") & (str_sub(this_colname, 16, 24) == "_GENEINFO")) {
    clinvar_date = as.character(str_sub(this_colname, 8, 15))
  }
}
if (clinvar_date != '') {
  #for (this_colname in colnames(data)) {
  for (i in 1 : length(colnames(data))) {
    this_colname = colnames(data)[i]
    clinvar_prefix = paste( 'CLINVAR', clinvar_date, '_', sep='' )
    if (str_sub(this_colname, 1, 16) == clinvar_prefix) {
      new_colname = paste( 'CLINVAR_', str_sub(this_colname, 17), sep='' )
      colnames(data)[i] = new_colname
    }
  }
}

# Rename the HGMD columns to remove date so that this and other programs know what the column names are, regardless of HGMD version.
hgmd_version = ''
for (this_colname in colnames(data)) {
  if ((str_sub(this_colname, 1, 7) == "HGMDPRO") & (str_sub(this_colname, 14) == "_GENEINFO")) {
    hgmd_version = as.character(str_sub(this_colname, 8, 13))
  } else {
    if ((str_sub(this_colname, 1, 7) == "HGMDPRO") & (str_sub(this_colname, 14) == "_GENE")) {
      hgmd_version = as.character(str_sub(this_colname, 8, 13))
    }
  }
}
if (hgmd_version != '') {
  for (i in 1 : length(colnames(data))) {
    this_colname = colnames(data)[i]
    hgmd_prefix = paste( 'HGMDPRO', hgmd_version, '_', sep='' )
    if (str_sub(this_colname, 1, 14) == hgmd_prefix) {
      new_colname = paste( 'HGMDPRO_', str_sub(this_colname, 15), sep='' )
      colnames(data)[i] = new_colname
    }
  }
}

# Rename sample_format fields to readable names
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

if (!('Allelic_depths_ref_alt' %in% colnames(data))) {
  data$Allelic_depths_ref_alt = "-1,-1"
  if (('FORMAT_NR' %in% colnames(data)) & ('FORMAT_NV' %in% colnames(data))) {
    data$Allelic_depths_ref_alt = paste( data$FORMAT_NR, data$FORMAT_NV, sep="," )
  }
}

# The following columns can contain gene names:
#	Gene.refGene=ABCB5				INFO_Gene.refGene
#	Gene.refGene==LINC02564\x3bLOC102723376
#	Gene.wgEncodeGencodeBasicV33=ABCB5		INFO_Gene.wgEncodeGencodeBasicV33
#	Gene.wgEncodeGencodeBasicV33=NONE\x3bAC215217.1
#	Gene.wgEncodeGencodeBasicV33=TUBB8\x3bZMYND11
#	GENE=ABCB5
#	HGMDPRO2019p4_GENE=ABCB5
#	CLINVAR20200629_GENEINFO=BRCA1:672
#	CLINVAR20200629_GENEINFO=CASQ2:845|VANGL1:81839
#	CLINVAR20200629_GENEINFO=TNNT2:7139
#	spliceai_gene=ABCB5
#	FANTOM5_CAGEr_humanHeart_transcriptStartSites=Name\x3dX237627570_237627644_+_240.4977261_LRRFIP1
#	FANTOM5_CAGEr_humanHeart_transcriptStartSites=Name\x3dX201365214_201365314_-_0_TNNT2,X201365186_201365281_+_202.7365015_AC119427.2

genes = read.table( genes_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )

# Extract gene name in columns Gene.refGene, Gene.wgEncodeGencodeBasicV33, and GENE.
# Multiple gene names are separated by \\x3b
data1 = data
data2 = data
if ("Gene.refGene" %in% colnames(data1)) {
  data2 = unique(as.data.frame( separate_rows(data1, Gene.refGene, sep="\\\\x3b", convert=FALSE) ))
  data1 = data2[(data2$Gene.refGene != "NONE"),]
}
if ("Gene.wgEncodeGencodeBasicV19" %in% colnames(data1)) {
  data2 = unique(as.data.frame( separate_rows(data1, Gene.wgEncodeGencodeBasicV19, sep="\\\\x3b", convert=FALSE) ))
  data1 = data2[(data2$Gene.wgEncodeGencodeBasicV19 != "NONE"),]
}
if ("Gene.wgEncodeGencodeBasicV33" %in% colnames(data1)) {
  data2 = unique(as.data.frame( separate_rows(data1, Gene.wgEncodeGencodeBasicV33, sep="\\\\x3b", convert=FALSE) ))
  data1 = data2[(data2$Gene.wgEncodeGencodeBasicV33 != "NONE"),]
}
if ("SPLICEOGEN2019_GENE" %in% colnames(data1)) {
  data2 = unique(as.data.frame( separate_rows(data1, SPLICEOGEN2019_GENE, sep="\\\\x3b", convert=FALSE) ))
  data1 = data2[(data2$SPLICEOGEN2019_GENE != "NONE"),]
}
if ("CLINVAR_GENEINFO" %in% colnames(data1)) {
  # At some point after here, when there is too much input, this program can crash due to running out of memory.
  # Change the split field delimiter of | to , because when separate_rows by |, it splits into a individual characters - probably a bug in separate_rows
  data1$CLINVAR_GENEINFO = str_replace_all(data1$CLINVAR_GENEINFO, fixed("|"), ",")
  data2 = unique(as.data.frame( separate_rows(data1, CLINVAR_GENEINFO, sep=",", convert=FALSE) ))
  data1 = data2[(data2$CLINVAR_GENEINFO != "NONE"),]
}
if ("FANTOM5_CAGEr_humanHeart_transcriptStartSites" %in% colnames(data1)) {
  data2 = unique(as.data.frame( separate_rows(data1, FANTOM5_CAGEr_humanHeart_transcriptStartSites, sep=",", convert=FALSE) ))
  data1 = data2[(data2$FANTOM5_CAGEr_humanHeart_transcriptStartSites != ""),]
  data2 = data1
}
rm(data1)

# Extract gene name in column CLINVAR_GENEINFO
# Gene name may be followed by :gene_number, eg. BRCA1:672
data3 = data2
data3$CLINVAR_GENEINFO_copy = data3$CLINVAR_GENEINFO
data3$CLINVAR_GENEINFO = ifelse( ( (grepl(":", data3$CLINVAR_GENEINFO, fixed=TRUE) == FALSE) & (data3$CLINVAR_GENEINFO != ".")), paste(data3$CLINVAR_GENEINFO, ":.", sep=""), data3$CLINVAR_GENEINFO )
data3$CLINVAR_GENEINFO = ifelse( (data3$CLINVAR_GENEINFO == ".") , ".:.", data3$CLINVAR_GENEINFO )
data3 = data3 %>% separate(CLINVAR_GENEINFO, c("CLINVAR_GENEINFO_gene_name", "CLINVAR_GENEINFO_gene_number"), sep=":", convert=FALSE)
rm(data2)

# Extract gene name in column FANTOM5_CAGEr_humanHeart_transcriptStartSites
# Gene name is the fifth component where components are separated by underscore
data4 = data3

if (nrow(data4) == 0) {
  write.table( data4, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
  quit()
}

if ("FANTOM5_CAGEr_humanHeart_transcriptStartSites" %in% colnames(data4)) {

  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = as.character(data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites)

  # The following are manipulation of data to extract gene name from FANTOM5 field so that it can be matched and the variant extracted by gene name.
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = ifelse( (data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites == '.'), "._._._._.", data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites )
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = gsub( "^Name\\\\x3dX", "", data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites )
  data4 = data4 %>% separate(FANTOM5_CAGEr_humanHeart_transcriptStartSites, c("FANTOM5_TSS_start", "FANTOM5_TSS_end", "FANTOM5_TSS_strand", "FANTOM5_TSS_tpm", "FANTOM5_TSS_gene"), sep="_")
}
rm(data3)

# extract variants that hit the genes of interest

got_gene = function(data, genes) {
  # data contains one row whose columns contain genes of a variant of a sample.
  # Most columns contain just one gene value. VEP_SYMBOL contains a bar-delimited list of genes.
  # genes contains the reference genes of interest.
  # This function works out whether any of the variant gene columns (only one row) contain a gene that is in the genes reference list.
  keep_this_row = FALSE
  for (i in 1:length(data)) {
    this_value = data[[i]]
    if (this_value %in% genes) {
      keep_this_row = TRUE
    }
    if (names(data)[i] == "VEP_SYMBOL") {
      if (grepl("|", this_value, fixed=TRUE)) {
        list_of_values = as.list(strsplit(this_value, '\\|')[[1]])
        for (this_vep_gene in list_of_values) {
          if (this_vep_gene %in% genes) {
            keep_this_row = TRUE
          }
        }
      }
    }
  }
  keep_this_row
}

list_gene_columns = c("Gene.refGene")
if ("Gene.wgEncodeGencodeBasicV19" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "Gene.wgEncodeGencodeBasicV19")
}
if ("Gene.wgEncodeGencodeBasicV33" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "Gene.wgEncodeGencodeBasicV33")
}
if ("sai_gene" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "sai_gene")
}
if ("saiI_gene" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "saiI_gene")
}
if ("spliceai_gene" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "spliceai_gene")
}
if ("spliceaiIndels_gene" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "spliceaiIndels_gene")
}
if ("spliceogen_gene" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "spliceogen_gene")
}
if ("CLINVAR_GENEINFO_gene_name" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "CLINVAR_GENEINFO_gene_name")
}
if ("GENE" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "GENE")
}
if ("SPLICEOGEN2019_GENE" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "SPLICEOGEN2019_GENE")
}
if ("HGMDPRO_GENE" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "HGMDPRO_GENE")
}
if ("FANTOM5_TSS_gene" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "FANTOM5_TSS_gene")
}
if ("PEXT_symbol" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "PEXT_symbol")
}
if ("VEP_SYMBOL" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "VEP_SYMBOL")
}
if ("SYMBOL" %in% colnames(data4)) {
  list_gene_columns = c(list_gene_columns, "SYMBOL")
}

data_gene_columns = data4[,list_gene_columns]
data5 = data4[ apply( data_gene_columns, 1, got_gene, genes$V1 ), ]
rm(data4)

if (nrow(data5) == 0) {
  write.table( data5, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
  quit()
}

data = data5
rm(data5)

data[is.na(data)] = ''
data[(data=='.')] = ''

# Remove any obviously superfluous fields.
#data$ANNOVAR_DATE = NULL

# Subset to variants that are possibly pathogenic. Fields to look at are:
#	Func.refGene
#	Func.wgEncodeGencodeBasicV33
#	HGMDPRO2019p4_CLASS
#	CLINVAR_CLNSIG
#	spliceai_site
#	FANTOM5_TSS_gene

# VCF FORMAT field of GATK SNPs usually contain GT:AD:DP:GQ:PL or GT:AD:DP:GQ:PGT:PID:PL
# When GATK3 ReadBackedPhasing is used to merge adjacent SNPs into MNP (multinucleotide polymorphism), FORMAT contains only GT:GQ
# This causes the separate command of the depths_ref_alt field to return NA. Set it to -1 to flag that it is lost depth information.
# Later, VPOT is parameterised to reject values < 10.

data$Allelic_depths_ref_alt_copy1 = data$Allelic_depths_ref_alt
data$Allelic_depths_ref_alt_copy2 = data$Allelic_depths_ref_alt
data$Allelic_depths_ref_alt_copy2 = ifelse( (data$Allelic_depths_ref_alt_copy2 == ""), "-1,-1", data$Allelic_depths_ref_alt_copy2 )
data1 = data %>% separate( Allelic_depths_ref_alt_copy2, c("Ref_depth", "Alt_depth"), sep=",", convert=FALSE)
data = data1
data$Allelic_depths_ref_alt = NULL
names(data)[names(data) == 'Allelic_depths_ref_alt_copy1'] = 'Allelic_depths_ref_alt'

# VEP consequences are listed at:
# https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
# In our data the VEP consequences for multiple transcripts appear in one field, separated by bar |
# VEP plugin UTRannotator consequences are listed at:
# https://github.com/ImperialCardioGenetics/UTRannotator

does_VEP_field_have_consequence_of_interest = function(data, list_of_deleterious_consequences) {

  VEP_consequence = unname(data[1]) # values are separated by bar |
  keep_this_row = FALSE
  VEP_consequence = gsub("\\&","\\|", VEP_consequence)
  data_consequence_list = unlist(strsplit(VEP_consequence,"\\|"))
  for (this_consequence in data_consequence_list) {
    if (this_consequence %in% list_of_deleterious_consequences) {
      keep_this_row = TRUE
      #print(paste(this_consequence, keep_this_row))
    }
  }
  #print(keep_this_row)
  keep_this_row
}

deleterious_VEP_consequences = c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion")
deleterious_UTRannotator_consequences = c("uAUG_gained", "uAUG_lost", "uSTOP_lost", "uSTOP gained", "uFrameShift")

# An example data$VEP_Consequence value:
# intron_variant&non_coding_transcript_variant|non_coding_transcript_exon_variant|intron_variant&non_coding_transcript_variant|intron_variant&non_coding_transcript_variant|non_coding_transcript_exon_variant|intron_variant&non_coding_transcript_variant

# An example of VEP_five_prime_UTR_variant_consequence:
# |uAUG_gained||||uAUG_gained|uAUG_gained||uAUG_gained||uAUG_gained||||||||||||||||uAUG_gained|uAUG_gained|uAUG_gained||||||||||||

# manta$gene_CDSexons_has_gene = as.numeric( lapply( manta$gene_CDSexons, FUN=genes_in_list, gene_list=gene_list ) )

data$extract_due_to_VEP_consequence = FALSE
if ("VEP_Consequence" %in% colnames(data)) {
  save_result = apply( data[,c("VEP_Consequence", "VEP_Consequence")], MARGIN=1, does_VEP_field_have_consequence_of_interest, deleterious_VEP_consequences)
  data$extract_due_to_VEP_consequence = ifelse( save_result, TRUE, data$extract_due_to_VEP_consequence )
}
if ("VEP_five_prime_UTR_variant_consequence" %in% colnames(data)) {
  save_result = apply(data[,c("VEP_five_prime_UTR_variant_consequence", "VEP_five_prime_UTR_variant_consequence")], MARGIN=1, does_VEP_field_have_consequence_of_interest, deleterious_UTRannotator_consequences)
  data$extract_due_to_VEP_consequence = ifelse( save_result, TRUE, data$extract_due_to_VEP_consequence )
}

data$extract_this_var = 0

if ("data$Func.refGene" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$Func.refGene %like% "exonic"), 1, data$extract_this_var )
  data$extract_this_var = ifelse( (data$Func.refGene %like% "splicing"), 1, data$extract_this_var )
}

if ("Func.wgEncodeGencodeBasicV19" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$Func.wgEncodeGencodeBasicV19 %like% "exonic"), 1, data$extract_this_var )
  data$extract_this_var = ifelse( (data$Func.wgEncodeGencodeBasicV19 %like% "splicing"), 1, data$extract_this_var )
}

if ("Func.wgEncodeGencodeBasicV33" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$Func.wgEncodeGencodeBasicV33 %like% "exonic"), 1, data$extract_this_var )
  data$extract_this_var = ifelse( (data$Func.wgEncodeGencodeBasicV33 %like% "splicing"), 1, data$extract_this_var )
}

if ("CLINVAR_CLNSIG" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$CLINVAR_CLNSIG %like% "athogenic"), 1, data$extract_this_var )
  data$extract_this_var = ifelse( (data$CLINVAR_CLNSIG %like% "ncertain"), 1, data$extract_this_var )
  data$extract_this_var = ifelse( (data$CLINVAR_CLNSIG %like% "risk"), 1, data$extract_this_var )
}

if ("HGMDPRO_CLASS" %in% colnames(data)) {
  data$extract_this_var = ifelse( !(data$HGMDPRO_CLASS %in% c("R",".","")), 1, data$extract_this_var )
}
if ("HGMDPRO_PROT" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$HGMDPRO_PROT %like% "p.M1"), 1, data$extract_this_var )
}

list_of_spliceai_columns = c("spliceai_acceptor_gain_delta_score", "spliceai_acceptor_loss_delta_score", "spliceai_donor_gain_delta_score", "spliceai_donor_loss_delta_score", "spliceaiIndels_acceptor_gain_delta_score", "spliceaiIndels_acceptor_loss_delta_score", "spliceaiIndels_donor_gain_delta_score", "spliceaiIndels_donor_loss_delta_score", "sai_AG_DS", "sai_AL_DS", "sai_DG_DS", "sai_DL_DS", "saiI_AG_DS", "saiI_AL_DS", "saiI_DG_DS", "saiI_DL_DS")

for (this_col in list_of_spliceai_columns) {
  if (this_col %in% colnames(data)) {
    data$temp_this_col = data[,c(this_col)]
    data$temp_this_col = ifelse( (data$temp_this_col %in% c(".","")), -1, data$temp_this_col )
    data$temp_this_col = as.numeric(as.character(data$temp_this_col))
    data$extract_this_var = ifelse( (data$temp_this_col >= 0.5), 1, data$extract_this_var )
    data$temp_this_col = NULL
  }
}

if ("FANTOM5_TSS_gene" %in% colnames(data)) {
  data$extract_this_var = ifelse( !(data$FANTOM5_TSS_gene %in% c(".","")), 1, data$extract_this_var )
}

if ("AAChange.refGene" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$AAChange.refGene %like% "p.M1"), 1, data$extract_this_var )
}
if ("AAChange.wgEncodeGencodeBasicV19" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$AAChange.wgEncodeGencodeBasicV19 %like% "p.M1"), 1, data$extract_this_var )
}
if ("AAChange.wgEncodeGencodeBasicV33" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$AAChange.wgEncodeGencodeBasicV33 %like% "p.M1"), 1, data$extract_this_var )
}

if ("mmsplice_delta_logit_psi" %in% colnames(data)) {
  data$extract_this_var = ifelse( (!(data$mmsplice_delta_logit_psi %in% c(".",""))), 1, data$extract_this_var )
}

if ("MutPredIndel_score" %in% colnames(data)) {
  data$extract_this_var = ifelse( !(data$MutPredIndel_score %in% c(".","")), 1, data$extract_this_var )
}

if ("extract_due_to_VEP_consequence" %in% colnames(data)) {
  data$extract_this_var = ifelse( (data$extract_due_to_VEP_consequence==TRUE), 1, data$extract_this_var )
}

patho = as.data.frame(data[ (data$extract_this_var==1) ,])

patho[is.na(patho)] = ''
patho$extract_due_to_VEP_consequence = NULL
patho$extract_this_var = NULL

write.table( patho, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


