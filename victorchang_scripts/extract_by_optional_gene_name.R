# module load R/3.6.1
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")

# Rscript extract_by_optional_gene_name.R infile outfile genes_file
# This program is similar to extract_by_gene_name.R,
# however extract_by_gene_name.R crashes when columns are not present
# whereas extract_by_optional_gene_name.R does not crash.

options(width=200)
options(stringsAsFactors = FALSE)
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( 'temp16.txt', 'tempout16.txt', 'list_genes_DCM_piezo1and2_withoutBars.txt' ) # for testing
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

  specific_input_is_empty = paste('Input file has no data: ',infile,sep='')
  input_is_empty = "CHROM\tPOS\tEND\tREF\tALT\tthere_are_no_variants"

  print( specific_input_is_empty )

  writeLines( input_is_empty, outfile )

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
#data2 = data
if ('Gene.refGene' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, Gene.refGene, sep="\\\\x3b", convert=FALSE) )
  data1 = data2[(data2$Gene.refGene != "NONE"),]
}
if ('Gene.wgEncodeGencodeBasicV33' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, Gene.wgEncodeGencodeBasicV33, sep="\\\\x3b", convert=FALSE) )
  data1 = data2[(data2$Gene.wgEncodeGencodeBasicV33 != "NONE"),]
}
if ('GENE' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, GENE, sep="\\\\x3b", convert=FALSE) )
  data1 = data2[(data2$GENE != "NONE"),]
}
if ('spliceogen_GENE' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, spliceogen_GENE, sep="\\\\x3b", convert=FALSE) )
  data1 = data2[(data2$spliceogen_GENE != "NONE"),]
}
if ('SPLICEOGEN2019_GENE' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, SPLICEOGEN2019_GENE, sep="\\\\x3b", convert=FALSE) )
  data1 = data2[(data2$SPLICEOGEN2019_GENE != "NONE"),]
}
if ('CLINVAR_GENEINFO' %in% colnames(data1)) {
  # At some point after here, when there is too much input, this program can crash due to running out of memory.
  # Change the split field delimiter of | to , because when separate_rows by |, it splits into a individual characters - probably a bug in separate_rows
  data1$CLINVAR_GENEINFO = str_replace_all(data1$CLINVAR_GENEINFO, fixed("|"), ",")
  data2 = unique( separate_rows(data1, CLINVAR_GENEINFO, sep=",", convert=FALSE) )
  data1 = data2[(data2$CLINVAR_GENEINFO != "NONE"),]
}
if ('FANTOM5_CAGEr_humanHeart_transcriptStartSites' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, FANTOM5_CAGEr_humanHeart_transcriptStartSites, sep=",", convert=FALSE) )
  data1 = data2[(data2$FANTOM5_CAGEr_humanHeart_transcriptStartSites != ""),]
  data2 = data1
}
if ('VEP_SYMBOL' %in% colnames(data1)) {
  data2 = unique( separate_rows(data1, VEP_SYMBOL, sep=",", convert=FALSE) )
  data1 = data2[(data2$VEP_SYMBOL != ""),]
  data2 = data1
}
rm(data1)

data3 = data2
if ('CLINVAR_GENEINFO' %in% colnames(data3)) {
  # Extract gene name in column CLINVAR_GENEINFO
  # Gene name may be followed by :gene_number, eg. BRCA1:672
  data3$CLINVAR_GENEINFO_copy = data3$CLINVAR_GENEINFO
  data3$CLINVAR_GENEINFO = ifelse( ( (grepl(":", data3$CLINVAR_GENEINFO, fixed=TRUE) == FALSE) & (data3$CLINVAR_GENEINFO != ".")), paste(data3$CLINVAR_GENEINFO, ":.", sep=""), data3$CLINVAR_GENEINFO )
  data3$CLINVAR_GENEINFO = ifelse( (data3$CLINVAR_GENEINFO == ".") , ".:.", data3$CLINVAR_GENEINFO )
  data3 = data3 %>% separate(CLINVAR_GENEINFO, c("CLINVAR_GENEINFO_gene_name", "CLINVAR_GENEINFO_gene_number"), sep=":", convert=FALSE)
}
rm(data2)

data4 = data3
if ('FANTOM5_CAGEr_humanHeart_transcriptStartSites' %in% colnames(data4)) {
  # Extract gene name in column FANTOM5_CAGEr_humanHeart_transcriptStartSites
  # Gene name is the fifth component where components are separated by underscore
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = as.character(data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites)

  # The following are manipulation of data to extract gene name from FANTOM5 field so that it can be matched and the variant extracted by gene name.
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = ifelse( (data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites == '.'), "._._._._.", data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites )
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = gsub( "^Name\\\\x3dX", "", data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites )
  data4 = data4 %>% separate(FANTOM5_CAGEr_humanHeart_transcriptStartSites, c("FANTOM5_TSS_start", "FANTOM5_TSS_end", "FANTOM5_TSS_strand", "FANTOM5_TSS_tpm", "FANTOM5_TSS_gene"), sep="_")
}
rm(data3)

# extract variants that hit the genes of interest

#got_gene = function(data, genes) {
#  gene1 = unname(data[1])
#  gene2 = data[2]
#  gene3 = data[3]
#  gene4 = data[4]
#  gene5 = data[5]
#  gene6 = data[6]
#  gene7 = data[7]
#  keep_this_row = FALSE
#  if ((gene1 %in% genes) | (gene2 %in% genes) | (gene3 %in% genes) | (gene4 %in% genes) | (gene5 %in% genes) | (gene6 %in% genes) | (gene7 %in% genes)) {
#    keep_this_row = TRUE
#  }
#  keep_this_row
#}
got_gene = function(data, genes) {
  gene1 = unname(data[1])
  keep_this_row = FALSE
  if (gene1 %in% genes) {
    keep_this_row = TRUE
  }
  if (names(data)[1] == "VEP_SYMBOL") {
    if (grepl("|", gene1, fixed=TRUE)) {
      list_of_values = as.list(strsplit(gene1, '\\|')[[1]])
      for (this_vep_gene in list_of_values) {
        if (this_vep_gene %in% genes) {
          keep_this_row = TRUE
        }
      }
    }
  }
  keep_this_row
}

data4$keep_this_row = 0
#data_gene_columns = data4[,c("Gene.refGene", "Gene.wgEncodeGencodeBasicV33", "GENE", "HGMDPRO2019p4_GENE", "spliceai_gene", "spliceaiIndels_gene", "spliceogen_gene", "CLINVAR_GENEINFO_gene_name", "FANTOM5_TSS_gene")]
#data5 = data4[ apply( data_gene_columns, 1, got_gene, genes$V1 ), ]
if ('Gene.refGene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("Gene.refGene", "Gene.refGene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('Gene.wgEncodeGencodeBasicV33' %in% colnames(data4)) {
  data_gene_columns = data4[,c("Gene.wgEncodeGencodeBasicV33", "Gene.wgEncodeGencodeBasicV33")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('GENE' %in% colnames(data4)) {
  data_gene_columns = data4[,c("GENE", "GENE")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('spliceogen_GENE' %in% colnames(data4)) {
  data_gene_columns = data4[,c("spliceogen_GENE", "spliceogen_GENE")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('SPLICEOGEN2019_GENE' %in% colnames(data4)) {
  data_gene_columns = data4[,c("SPLICEOGEN2019_GENE", "SPLICEOGEN2019_GENE")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('HGMDPRO2019p4_GENE' %in% colnames(data4)) {
  data_gene_columns = data4[,c("HGMDPRO2019p4_GENE", "HGMDPRO2019p4_GENE")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('HGMDPRO2020p3_GENE' %in% colnames(data4)) {
  data_gene_columns = data4[,c("HGMDPRO2020p3_GENE", "HGMDPRO2020p3_GENE")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('sai_gene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("sai_gene", "sai_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('saiI_gene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("saiI_gene", "saiI_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('spliceai_gene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("spliceai_gene", "spliceai_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('spliceaiIndels_gene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("spliceaiIndels_gene", "spliceaiIndels_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('spliceogen_gene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("spliceogen_gene", "spliceogen_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('CLINVAR_GENEINFO_gene_name' %in% colnames(data4)) {
  data_gene_columns = data4[,c("CLINVAR_GENEINFO_gene_name", "CLINVAR_GENEINFO_gene_name")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('FANTOM5_TSS_gene' %in% colnames(data4)) {
  data_gene_columns = data4[,c("FANTOM5_TSS_gene", "FANTOM5_TSS_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('FANTOM5_CAGEr_humanHeart_transcriptStartSites' %in% colnames(data4)) {
  # Extract gene name in column FANTOM5_CAGEr_humanHeart_transcriptStartSites
  # Gene name is the fifth component where components are separated by underscore
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = as.character(data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites)
 # The following are manipulation of data to extract gene name from FANTOM5 field so that it can be matched and the variant extracted by gene name.
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = ifelse( (data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites == '.'), "._._._._.", data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites )
  data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites = gsub( "^Name\\\\x3dX", "", data4$FANTOM5_CAGEr_humanHeart_transcriptStartSites )
  data4 = data4 %>% separate(FANTOM5_CAGEr_humanHeart_transcriptStartSites, c("FANTOM5_TSS_start", "FANTOM5_TSS_end", "FANTOM5_TSS_strand", "FANTOM5_TSS_tpm", "FANTOM5_TSS_gene"), sep="_")
  data_gene_columns = data4[,c("FANTOM5_TSS_gene", "FANTOM5_TSS_gene")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}
if ('VEP_SYMBOL' %in% colnames(data4)) {
  data_gene_columns = data4[,c("VEP_SYMBOL", "VEP_SYMBOL")]
  data4$keep_this_row = ifelse( apply( data_gene_columns, 1, got_gene, genes$V1 ), 1, data4$keep_this_row )
}

data5 = data4[ (data4$keep_this_row==1), ]
data5$keep_this_row = NULL
rm(data4)

if (nrow(data5) == 0) {
  write.table( data5, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
} else {

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

if ('Allelic_depths_ref_alt' %in% colnames(data)) {
  data$Allelic_depths_ref_alt_copy1 = data$Allelic_depths_ref_alt
  data$Allelic_depths_ref_alt_copy2 = data$Allelic_depths_ref_alt
  data1 = data %>% separate( Allelic_depths_ref_alt_copy2, c("Ref_depth", "Alt_depth"), sep=",", convert=FALSE)
  data = data1
  data$Allelic_depths_ref_alt = NULL
  names(data)[names(data) == 'Allelic_depths_ref_alt_copy1'] = 'Allelic_depths_ref_alt'
}

patho = data # This script was copied from extract_pathogenic_by_gene_name and is the same except for this line.
patho[is.na(patho)] = ''

write.table( patho, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )

}

