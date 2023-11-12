# module load R/3.6.1
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")
# export R_LIBS_USER=/g/data/jb96/software/victorchang_scripts/R_libraries_for_victorchang_scripts

# Rscript add_gene_flag_to_report.R infile outfile gene_flag gene_names_file gene_regions_file report_gene_column_names

#Sys.getenv('R_LIBS_USER')
#.libPaths()

options(width=200)
options(stringsAsFactors = FALSE)
options(scipen=10000)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'A020029517111.cnvnator.filter.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.bndInsideExon.compareSamples.extractTGAgenes.tsv', 'A020029517111.cnvnator.flag.tsv', 'NDD', 'list_genes_tga_NDD_659_genes_withoutBars.txt', 'list_genes_tga_NDD_659_genes_RefSeq_regions_tab_delimited.txt', 'gene_entireGene,gene_CDSexons' ) # for testing

infile = as.character(args[1])
outfile = as.character(args[2])
gene_flag = as.character(args[3])
gene_names_file = as.character(args[4])
gene_regions_file = as.character(args[5])
report_gene_column_names = as.character(args[6])

#print('report_gene_column_names')
#print(report_gene_column_names)

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

colnames(data)[which(colnames(data) == 'X.CHROM')] = 'chrom'
colnames(data)[which(colnames(data) == '#CHROM')] = 'chrom'
colnames(data)[which(colnames(data) == 'CHROM')] = 'chrom'
colnames(data)[which(colnames(data) == 'X.Chr')] = 'chrom'
colnames(data)[which(colnames(data) == '#Chr')] = 'chrom'
colnames(data)[which(colnames(data) == 'Chr')] = 'chrom'
colnames(data)[which(colnames(data) == 'X.Chrom')] = 'chrom'
colnames(data)[which(colnames(data) == '#Chrom')] = 'chrom'
colnames(data)[which(colnames(data) == 'Chrom')] = 'chrom'
colnames(data)[which(colnames(data) == 'X.chrom')] = 'chrom'
colnames(data)[which(colnames(data) == '#chrom')] = 'chrom'

#print('colnames(data)')
#print(colnames(data))

genes1 = read.table( gene_names_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(genes1) = c("gene_for_matching")
genes = unique( paste( "|", genes1$gene, "|", sep="" ) )

if ((gene_regions_file != '.') & (gene_regions_file != '')) {

  regions = read.table( gene_regions_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
  colnames(regions) = c( "chrom", "start", "end" )
  regions$chrom = as.character(regions$chrom)
  regions$start = as.numeric(regions$start)
  regions$end = as.numeric(regions$end)
}

match_by_gene_names = function(data_row, genes_list, output) {
  match_element = data_row[1]
  found = FALSE
  for (this_gene in genes_list) {
    if (grepl(this_gene, match_element, fixed=TRUE)) {
      found = TRUE
    }
  }
  found
}

match_by_gene_regions = function(data_row, regions_table, output) {
  match_chrom = as.character(data_row[1])
  match_start = as.numeric(data_row[2])
  match_end = as.numeric(data_row[3])
  found = FALSE
  regions_table_for_one_chrom = regions_table[ (regions_table$chrom==match_chrom), ]
  table_of_matches = regions_table_for_one_chrom[ ((regions_table_for_one_chrom$start<=match_start)&(regions_table_for_one_chrom$end>=match_start)) | ((regions_table_for_one_chrom$start<=match_end)&(regions_table_for_one_chrom$end>=match_end)) | ((regions_table_for_one_chrom$start<=match_start)&(regions_table_for_one_chrom$end>=match_end)) | ((regions_table_for_one_chrom$start>=match_start)&(regions_table_for_one_chrom$end<=match_end)), ]
  if (nrow(table_of_matches) > 0) {
    found = TRUE
  }
  found
}

data2 = data
data2$new_column_containing_flag = ""

list_report_gene_column_names = str_split( report_gene_column_names, "," ) # [1] "gene_entireGene" "gene_CDSexons"

#print('colnames(data2)')
#print(colnames(data2))
#print('list_report_gene_column_names')
#print(list_report_gene_column_names)

for (this_report_column_name in list_report_gene_column_names[[1]]) {

  if (this_report_column_name %in% colnames(data2)) {

    colnum = match( this_report_column_name, colnames(data2) )

    data2$gene_for_matching = data2[ ,colnum ]
    data2$gene_for_matching = paste( "|", data2$gene_for_matching, "|", sep="" )
    data2$gene_for_matching = gsub( ", ", "|, |", data2$gene_for_matching )

    data3 = data2[, c(ncol(data2),1:(ncol(data2)-1)) ]
    # mictro - 9/10/23 ; commented out EI2023 change
    # adding this inside the for loop will overwrite already assigned "YES" for previous match
    # This can cause a "." to be in the column instead of "YES"
    #data3$new_column_containing_flag = "." #EI2023 
    data3$new_column_containing_flag = ifelse( apply( data3, 1, match_by_gene_names, genes ), "YES", data3$new_column_containing_flag )

    data3$gene_for_matching = NULL
    data2 = data3
  }
}

data3 = unique(data2)

col1 = 1
col2 = 2
col3 = 3
if ('chrom' %in% colnames(data3)) {
  true_chrom_col = which(colnames(data3) == 'chrom')
  col1 = true_chrom_col
  col2 = col1 + 1
  col3 = col2 + 1
}
if (colnames(data3)[col3] == 'ID') {
  col3 = col2
}

# Put a key of chrom,pos,end at the front of dataframe, for use with regions filtering.
key = data3[,c(col1,col2,col3)]
colnames(key) = c( "temp_key_chrom", "temp_key_start", "temp_key_end" )
new_data3 = cbind( key, data3 )
data3 = new_data3

#EI2023 data3$new_column_containing_flag_2 = ""
data3$new_column_containing_flag_2 = "." #EI2023 
data3[,1] = as.character(data3[,1])
data3[,2] = as.numeric(data3[,2])
data3[,3] = as.numeric(data3[,3])

if ((gene_regions_file != '.') & (gene_regions_file != '')) {

  data3$new_column_containing_flag_2 = ifelse( apply( data3, 1, match_by_gene_regions, regions ), "YES", data3$new_column_containing_flag_2 )
}

data3$new_column_containing_flag = ifelse( (data3$new_column_containing_flag_2=="YES"), "YES", data3$new_column_containing_flag )
# mictro 9/10/23 - default to "." for new_column_containing_flag - if value is ""
data3$new_column_containing_flag = ifelse( (data3$new_column_containing_flag==""), ".", data3$new_column_containing_flag )
data3$new_column_containing_flag_2 = NULL
data3$temp_key_chrom = NULL
data3$temp_key_start = NULL
data3$temp_key_end = NULL

names(data3)[names(data3) == 'new_column_containing_flag'] = gene_flag

data3[is.na(data3)] = ""

write.table( data3, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


