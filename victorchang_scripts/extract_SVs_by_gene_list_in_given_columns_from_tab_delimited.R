# module load R/3.6.1
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")
# export R_LIBS_USER=/g/data/jb96/software/victorchang_scripts/R_libraries_for_victorchang_scripts

# Rscript extract_SVs_by_gene_list_in_given_columns_from_tab_delimited.R 19W001100.manta.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.collapsed.tsv 19W001100.manta.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.collapsed.CHD_tier1_106_genes_plus_extended_505_genes.tsv list_genes_CHD_tier1_106_genes_plus_extended_505_genes_withoutBars.txt

#Sys.getenv('R_LIBS_USER')
#.libPaths()

options(width=200)
options(stringsAsFactors = FALSE)
library(data.table)
library(dplyr)
library(tidyr)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'sample.gridss.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.tsv',  'sample.gridss.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.DCM_tier1.tsv', 'list_genes_DCM_tier1_list_of_32_genes_withoutBars.txt', 'gene_CDSexons,gene_entireGene,HGMDPRO2019p4_GENE,FANTOM5_heart_TSS_start_end_strand_tpm_gene' ) # for testing

infile = as.character(args[1])
outfile = as.character(args[2])
genefile = as.character(args[3])
report_gene_column_names = as.character(args[4])

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

  specific_input_is_empty = paste('Input file has header but no data: ',infile,sep='')
  input_is_empty = "CHROM\tPOS\tEND\tREF\tALT\tthere_are_no_variants"

  print( specific_input_is_empty )

  writeLines( input_is_empty, outfile )

  quit()
}

genes1 = read.table( genefile, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(genes1) = c("gene")
genes = unique( paste( "|", genes1$gene, "|", sep="" ) )

# we will search for the genes that appear in the following columns
gene_columns = strsplit(report_gene_column_names, ",", fixed=TRUE)[[1]]

# subset data to those having genes present, so that R is working on a smaller amount of data
data$temp_has_gene = ""
for (this_gene_col_name in gene_columns) {
  if (this_gene_col_name %in% colnames(data)) {
    this_gene_col_num = which( colnames(data)==this_gene_col_name )
    data$temp_has_gene = ifelse( (data[,this_gene_col_num] != "."), "yes", data$temp_has_gene )
  }
}
data2 = data[(data$temp_has_gene=="yes"),]
data2$temp_has_gene = NULL

extract_fantom5_heart_tss_genes = function( data_row, fantom5_col_num ) {
  out_data = ''
  fantom5_data = data_row[fantom5_col_num]
  arr1 = strsplit(fantom5_data, split = ", ")
  arr2 = arr1[[1]]
  for (this_element in arr2) {
    arr3 = strsplit(this_element, split = "_")
    arr4 = arr3[[1]]
    if (!is.na(arr4[5])) {
      this_gene = arr4[5]
      if (this_gene != '') {
        if (out_data == '') {
          out_data = this_gene
        } else {
          new_out_data = paste( out_data, this_gene, sep=", " )
          out_data = new_out_data
        }
      }
    }
  }
  out_data
}

# Column FANTOM5_heart_TSS_start_end_strand_tpm_gene has gene names embedded as the 5th underscore-delimited column,
# and there can be multiple of these 5-column values separated by comma-and-space.
# Extract the genes as comma-and-space-delimited so that they can be compared to the genes of interest, 
# put into new column called emp_fantom5_heart_tss.
if (("FANTOM5_heart_TSS_start_end_strand_tpm_gene" %in% gene_columns) & ("FANTOM5_heart_TSS_start_end_strand_tpm_gene" %in% colnames(data2))) {
  fantom5_col_num = which( colnames(data2)=="FANTOM5_heart_TSS_start_end_strand_tpm_gene" )
  data2$temp_fantom5_heart_tss = apply( data2, 1, extract_fantom5_heart_tss_genes, fantom5_col_num )
}

# Now compare the genes in the data to the genes of interest.
data2$temp_has_gene = ""
for (this_gene_col_name in gene_columns) {
  if (this_gene_col_name != "FANTOM5_heart_TSS_start_end_strand_tpm_gene") { # will use the new column temp_fantom5_heart_tss instead
    this_gene_col_num = which( colnames(data)==this_gene_col_name )
    data2$match1 = paste( "|", data2[,this_gene_col_num], "|", sep="" )
    data2$match1 = gsub( ", ", "|, |", data2$match1 )
    for (this_gene in genes) {
      data2$temp_has_gene = ifelse( (grepl(this_gene, data2$match1, fixed=TRUE)), "yes", data2$temp_has_gene )
    }
  }
}
data2$temp_fantom5_heart_tss = NULL

data3 = data2[(data2$temp_has_gene=="yes"),]
data3$temp_has_gene = NULL
data3$match1 = NULL

data4 = unique(data3)
data4[is.na(data4)] = ""

write.table( data4, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


