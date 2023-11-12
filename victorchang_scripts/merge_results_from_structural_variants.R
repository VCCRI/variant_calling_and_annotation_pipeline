
# This program assumes that there is at least one matching variant between manta and gridss.
# It then looks to see if cnvnator variants match any of those manta_matches_gridss variants.
# This program outputs the manta_matches_gridss variants, and any cnvnator variants that match them.

options(width=200)
options(stringsAsFactors = FALSE)
library(tidyr)
library(dplyr)
library(data.table)
library(sqldf)

args = commandArgs(trailingOnly=TRUE) # for production

#args=c( '19W001092.manta_gridss_cnvnator.merge.CHD_tier1_106_genes_plus_extended_505_genes_plus_CFAP54_CFC1', '19W001092.manta.CHD_tier1_106_genes_plus_extended_505_genes_plus_CFAP54_CFC1.tsv', '19W001092.gridss_1_4_1.CHD_tier1_106_genes_plus_extended_505_genes_plus_CFAP54_CFC1.tsv', '19W001092.cnvnator.CHD_tier1_106_genes_plus_extended_505_genes_plus_CFAP54_CFC1.tsv', 'list_genes_CHD_tier1_list_of_106_genes_withoutBars.txt', 'CHD_tier1_list_of_106_genes' )
#args=c( '19W001300.manta_gridss_cnvnator.merge.DCM_tier1_list_of_32_genes', '19W001300.manta.DCM_tier1_list_of_32_genes.tsv', '19W001300.gridss_1_4_1.DCM_tier1_list_of_32_genes.tsv', '19W001300.cnvnator.DCM_tier1_list_of_32_genes.tsv', '.', '.' )

#args=c( '19F00552.manta_gridss_cnvnator.merge', '19F00552.gridss_1_4_1.DCM_allJun2020_585_genes.vpot_final_output_file.vpot_gt_0.txt', '19F00552.manta.DCM_allJun2020_585_genes.vpot_final_output_file.vpot_gt_0.txt', '19F00552.cnvnator.DCM_allJun2020_585_genes.vpot_final_output_file.vpot_gt_0.txt', '.', '.' )

#args=c( '/g/data/a32/quarterly_x1_10TB/WGS/BAM/AGHA/AGHA_4/working_directory/merge_results/', '/g/data/a32/quarterly_x1_10TB/WGS/BAM/AGHA/AGHA_4/working_directory/manta_annotate/19F00552.manta.DCM_allJun2020_585_genes.vpot_final_output_file.vpot_gt_0.txt', '/g/data/a32/quarterly_x1_10TB/WGS/BAM/AGHA/AGHA_4/working_directory/gridss_annotate/19F00552.gridss_1_4_1.DCM_allJun2020_585_genes.vpot_final_output_file.vpot_gt_0.txt', '/g/data/a32/quarterly_x1_10TB/WGS/BAM/AGHA/AGHA_4/working_directory/cnvnator_annotate/19F00552.cnvnator.DCM_allJun2020_585_genes.vpot_final_output_file.vpot_gt_0.txt', '.', '.' )

#args=c( '19W001402.manta_gridss_cnvnator.merge.CHD_allJun2020_888_genes.tsv', '19W001402.manta.CHD_allJun2020_888_genes.tsv', '19W001402.gridss_1_4_1.CHD_allJun2020_888_genes.tsv', '19W001402.cnvnator.CHD_allJun2020_888_genes.tsv', '.', '.' )

#args=c("19F00553.manta_gridss_cnvnator.merge.DCM_allJun2020_585_genes", "19F00553.manta.DCM_allJun2020_585_genes.tsv", "19F00553.gridss_1_4_1.DCM_allJun2020_585_genes.tsv", "19F00553.cnvnator.DCM_allJun2020_585_genes.tsv", ".", ".")

outfile_prefix = as.character(args[1])
manta_file = as.character(args[2])
gridss_file = as.character(args[3])
cnvnator_file = as.character(args[4])
gene_list_file = as.character(args[5])
gene_list_id = as.character(args[6])

outfile_all = paste( outfile_prefix, ".all.tsv", sep="" )
outfile_truepos = paste( outfile_prefix, ".truePos.tsv", sep="" )
outfile_truepos_exons = paste( outfile_prefix, ".truePos_exons.tsv", sep="" )
outfile_truepos_splicing = paste( outfile_prefix, ".truePos_splicing.tsv", sep="" )

manta = read.table( manta_file, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
gridss = read.table( gridss_file, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
if (cnvnator_file != '.') {
  cnvnator = read.table( cnvnator_file, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
}

names(manta)[names(manta) == 'CHROM'] = 'chrom'
names(manta)[names(manta) == 'X.CHROM'] = 'chrom'
names(manta)[names(manta) == 'START'] = 'start'
names(manta)[names(manta) == 'POS'] = 'start'
names(manta)[names(manta) == 'END'] = 'end'
names(manta)[names(manta) == 'REF'] = 'ref'
names(manta)[names(manta) == 'Ref'] = 'ref'
names(manta)[names(manta) == 'ALT'] = 'alt'
names(manta)[names(manta) == 'Alt'] = 'alt'
names(manta)[names(manta) == 'SVLEN'] = 'svlen'
names(manta)[names(manta) == 'SVTYPE'] = 'svtype'
names(manta)[names(manta) == 'gnomadsv_SVTYPE_AF'] = 'gnomadsv_svtype_af'
manta$chrom = as.character(manta$chrom)
manta$start = as.numeric(manta$start)
manta$end = as.numeric(manta$end)
manta$svlen = as.numeric(manta$svlen)

names(gridss)[names(gridss) == 'CHROM'] = 'chrom'
names(gridss)[names(gridss) == 'X.CHROM'] = 'chrom'
names(gridss)[names(gridss) == 'START'] = 'start'
names(gridss)[names(gridss) == 'POS'] = 'start'
names(gridss)[names(gridss) == 'END'] = 'end'
names(gridss)[names(gridss) == 'REF'] = 'ref'
names(gridss)[names(gridss) == 'Ref'] = 'ref'
names(gridss)[names(gridss) == 'ALT'] = 'alt'
names(gridss)[names(gridss) == 'Alt'] = 'alt'
names(gridss)[names(gridss) == 'SVLEN'] = 'svlen'
names(gridss)[names(gridss) == 'SVTYPE'] = 'svtype'
names(gridss)[names(gridss) == 'gnomadsv_SVTYPE_AF'] = 'gnomadsv_svtype_af'
gridss$chrom = as.character(gridss$chrom)
gridss$start = as.numeric(gridss$start)
gridss$end = as.numeric(gridss$end)
gridss$svlen = as.numeric(gridss$svlen)

we_have_cnvnator = "no"
we_have_manta_cnvnator = "no"
we_have_gridss_cnvnator = "no"

if (cnvnator_file != '.') {
  if (nrow(cnvnator) > 0) {
    we_have_cnvnator = "yes"
    names(cnvnator)[names(cnvnator) == 'CHROM'] = 'chrom'
    names(cnvnator)[names(cnvnator) == 'X.CHROM'] = 'chrom'
    names(cnvnator)[names(cnvnator) == 'POS'] = 'start'
    names(cnvnator)[names(cnvnator) == 'START'] = 'start'
    names(cnvnator)[names(cnvnator) == 'END'] = 'end'
    names(cnvnator)[names(cnvnator) == 'REF'] = 'ref'
    names(cnvnator)[names(cnvnator) == 'Ref'] = 'ref'
    names(cnvnator)[names(cnvnator) == 'ALT'] = 'alt'
    names(cnvnator)[names(cnvnator) == 'Alt'] = 'alt'
    names(cnvnator)[names(cnvnator) == 'SVLEN'] = 'svlen'
    names(cnvnator)[names(cnvnator) == 'SVTYPE'] = 'svtype'
    names(cnvnator)[names(cnvnator) == 'gnomadsv_SVTYPE_AF'] = 'gnomadsv_svtype_af'
    cnvnator$chrom = as.character(cnvnator$chrom)
    cnvnator$start = as.numeric(cnvnator$start)
    cnvnator$end = as.numeric(cnvnator$end)
  }
}

genes_in_list <- function( params, gene_list ) {
  list_genes = params[1]
  list2 = strsplit( list_genes, ", " )
  genes_in_sample = unlist(list2)
  result = FALSE
  for( this_gene in gene_list ) {
    if (this_gene %in% genes_in_sample) {
      result = TRUE
    }
  }
  result
}

if (gene_list_file != ".") {

  gene_list_table = read.table( gene_list_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
  colnames(gene_list_table) = c("gene")
  gene_list = gene_list_table$gene

  manta$gene_CDSexons_has_gene = as.numeric( lapply( manta$gene_CDSexons, FUN=genes_in_list, gene_list=gene_list ) )
  manta$gene_entireGene_has_gene = as.numeric( lapply( manta$gene_entireGene, FUN=genes_in_list, gene_list=gene_list ) )
  colnames(manta)[colnames(manta) == "gene_CDSexons_has_gene"] = paste( "gene_CDSexons_has_", gene_list_id, sep="" )
  colnames(manta)[colnames(manta) == "gene_entireGene_has_gene"] = paste( "gene_entireGene_has_", gene_list_id, sep="" )

  gridss$gene_CDSexons_has_gene = as.numeric( lapply( gridss$gene_CDSexons, FUN=genes_in_list, gene_list=gene_list ) )
  gridss$gene_entireGene_has_gene = as.numeric( lapply( gridss$gene_entireGene, FUN=genes_in_list, gene_list=gene_list ) )
  colnames(gridss)[colnames(gridss) == "gene_CDSexons_has_gene"] = paste( "gene_CDSexons_has_", gene_list_id, sep="" )
  colnames(gridss)[colnames(gridss) == "gene_entireGene_has_gene"] = paste( "gene_entireGene_has_", gene_list_id, sep="" )

  if (we_have_cnvnator == "yes") {
    cnvnator$gene_CDSexons_has_gene = as.numeric( lapply( cnvnator$gene_CDSexons, FUN=genes_in_list, gene_list=gene_list ) )
    cnvnator$gene_entireGene_has_gene = as.numeric( lapply( cnvnator$gene_entireGene, FUN=genes_in_list, gene_list=gene_list ) )
    colnames(cnvnator)[colnames(cnvnator) == "gene_CDSexons_has_gene"] = paste( "gene_CDSexons_has_", gene_list_id, sep="" )
    colnames(cnvnator)[colnames(cnvnator) == "gene_entireGene_has_gene"] = paste( "gene_entireGene_has_", gene_list_id, sep="" )
  }
}

manta1 = manta[(manta$svtype != "BND"),] 
manta = manta1

manta$abs_svlen = as.numeric(as.character(manta$svlen))
manta$abs_svlen = abs(manta$svlen)
manta$abs_svlen = ifelse( is.na(manta$abs_svlen), 0, manta$abs_svlen )
manta$abs_svlen = ifelse( (manta$abs_svlen==0), abs(manta$end - manta$start), manta$abs_svlen )
manta$abs_svlen_15percent = round( manta$abs_svlen * 0.15, digits=0 )
manta$abs_svlen_20percent = round( manta$abs_svlen * 0.20, digits=0 )
manta$abs_svlen_70percent = round( manta$abs_svlen * 0.70, digits=0 )
colnames(manta) = paste( "manta_", colnames(manta), sep="" )

gridss$abs_svlen = as.numeric(as.character(gridss$svlen))
gridss$abs_svlen = abs(gridss$svlen)
gridss$abs_svlen = ifelse( is.na(gridss$abs_svlen), 0, gridss$abs_svlen )
gridss$abs_svlen = ifelse( (gridss$abs_svlen==0), abs(gridss$end - gridss$start), gridss$abs_svlen )
gridss$abs_svlen_15percent = round( gridss$abs_svlen * 0.15, digits=0 )
gridss$abs_svlen_20percent = round( gridss$abs_svlen * 0.20, digits=0 )
gridss$abs_svlen_70percent = round( gridss$abs_svlen * 0.70, digits=0 )
colnames(gridss) = paste( "gridss_", colnames(gridss), sep="" )

if (we_have_cnvnator == "yes") {
  cnvnator$abs_svlen = as.numeric(as.character(cnvnator$svlen))
  cnvnator$abs_svlen = cnvnator$end - cnvnator$start + 1
  cnvnator$abs_svlen = ifelse( is.na(cnvnator$abs_svlen), 0, cnvnator$abs_svlen )
  cnvnator$abs_svlen = ifelse( (cnvnator$abs_svlen==0), abs(cnvnator$end - cnvnator$start), cnvnator$abs_svlen )
  cnvnator$abs_svlen_15percent = round( cnvnator$abs_svlen * 0.15, digits=0 )
  cnvnator$abs_svlen_20percent = round( cnvnator$abs_svlen * 0.20, digits=0 )
  cnvnator$abs_svlen_70percent = round( cnvnator$abs_svlen * 0.70, digits=0 )
  colnames(cnvnator) = paste( "cnvnator_", colnames(cnvnator), sep="" )
}

# To be a matching variant, need to have 70% overlap, and each end to be within 20% of each other.

manta1 = as.data.frame(manta[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "manta_abs_svlen_15percent", "manta_abs_svlen_20percent", "manta_abs_svlen_70percent")])
gridss1 = as.data.frame(gridss[,c("gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt", "gridss_abs_svlen_15percent", "gridss_abs_svlen_20percent", "gridss_abs_svlen_70percent")])

merge1 = sqldf('select * from manta1 m, gridss1 g where (m.manta_chrom=g.gridss_chrom) and (abs(m.manta_start-g.gridss_start)<=min(m.manta_abs_svlen_20percent,g.gridss_abs_svlen_20percent)) and (abs(m.manta_end-g.gridss_end)<=min(m.manta_abs_svlen_20percent,g.gridss_abs_svlen_20percent)) and ((min(m.manta_end,g.gridss_end)-min(m.manta_start,g.gridss_start))>=max(m.manta_abs_svlen_70percent,g.gridss_abs_svlen_70percent))')
# merge1 contains all manta variants on a chromosome paired with all gridss variants on the same chromosome, and some of their data columns

merge2 = unique( merge1[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt")] )
# merge2 contains variables in both manta and gridss, and only their key columns

if (we_have_cnvnator == "yes") {

  cnvnator3 = cnvnator[,c("cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size", "cnvnator_abs_svlen_15percent", "cnvnator_abs_svlen_20percent", "cnvnator_abs_svlen_70percent")]

  merge3 = sqldf('select * from merge1 m, cnvnator3 c where (m.manta_chrom=c.cnvnator_chrom) and (   ((abs(m.manta_start-c.cnvnator_start)<=min(m.manta_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and (abs(m.manta_end-c.cnvnator_end)<=min(m.manta_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and ((min(m.manta_end,c.cnvnator_end)-min(m.manta_start,c.cnvnator_start))>=max(m.manta_abs_svlen_70percent,c.cnvnator_abs_svlen_70percent)))   or   ((abs(m.gridss_start-c.cnvnator_start)<=min(m.gridss_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and (abs(m.gridss_end-c.cnvnator_end)<=min(m.gridss_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and ((min(m.gridss_end,c.cnvnator_end)-min(m.gridss_start,c.cnvnator_start))>=max(m.gridss_abs_svlen_70percent,c.cnvnator_abs_svlen_70percent)))   )')
  # merge3 contains variables in all 3 of manta and gridss and cnvnator, and some of their data columns

  merge4 = unique( merge3[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt", "cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size")] )
  # merge4 contains variables in all 3 of manta and gridss and cnvnator, and only their key columns

  merge5 = merge( x=merge4, y=merge2, by=c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt"), all.x=TRUE, all.y=TRUE )
  # merge5 contains variables in all 3 of manta and gridss and cnvnator, and variables in only manta and gridss, and only the key columns of all 3 data sources

  merge5$cnvnator_chrom = ifelse( is.na(merge5$cnvnator_chrom), '', merge5$cnvnator_chrom )
  merge5$cnvnator_start = ifelse( is.na(merge5$cnvnator_start), 0, merge5$cnvnator_start )
  merge5$cnvnator_end = ifelse( is.na(merge5$cnvnator_end), 0, merge5$cnvnator_end )
  merge5$cnvnator_cnv_type = ifelse( is.na(merge5$cnvnator_cnv_type), '', merge5$cnvnator_cnv_type )
  merge5$cnvnator_cnv_size = ifelse( is.na(merge5$cnvnator_cnv_size), 0, merge5$cnvnator_cnv_size )
}

if (we_have_cnvnator == "no") {
  merge5_1 = merge2

  merge5_1$cnvnator_chrom = ''
  merge5_1$cnvnator_start = 0
  merge5_1$cnvnator_end = 0
  merge5_1$cnvnator_cnv_type = ''
  merge5_1$cnvnator_cnv_size = 0

  merge5_2 = merge5_1[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt", "cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size")]
}

if (we_have_cnvnator == "yes") {

  merge6_1 = sqldf('select * from manta1 m, cnvnator3 c where (m.manta_chrom=c.cnvnator_chrom) and (   ((abs(m.manta_start-c.cnvnator_start)<=min(m.manta_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and (abs(m.manta_end-c.cnvnator_end)<=min(m.manta_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and ((min(m.manta_end,c.cnvnator_end)-min(m.manta_start,c.cnvnator_start))>=max(m.manta_abs_svlen_70percent,c.cnvnator_abs_svlen_70percent)))   )')

  # merge6_1 contains variables in manta and cnvnator, regardless of whether they exist in gridss, and some of their data columns

  merge6_2 = unique( merge6_1[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size")] )
  # merge6_2 contains variables in manta and cnvnator, regardless of whether they exist in gridss, and only their key columns

  merge6_3 = sqldf('select * from merge6_2 t1 where not exists (select * from merge2 t2 where (t1.manta_chrom=t2.manta_chrom and t1.manta_start=t2.manta_start and t1.manta_end=t2.manta_end and t1.manta_ref=t2.manta_ref and t1.manta_alt=t2.manta_alt))')
  # merge6_3 contains variables in manta and cnvnator that do not exist in gridss, and only their key columns

  if (nrow(merge6_3) > 0) {

    we_have_manta_cnvnator = "yes"

    merge6_3$gridss_chrom = ''
    merge6_3$gridss_start = 0
    merge6_3$gridss_end = 0
    merge6_3$gridss_ref = ''
    merge6_3$gridss_alt = ''

    merge6_4 = merge6_3[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt", "cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size")]
  }
}

if (we_have_cnvnator == "yes") {

  merge7_1 = sqldf('select * from gridss1 m, cnvnator3 c where (m.gridss_chrom=c.cnvnator_chrom) and (   ((abs(m.gridss_start-c.cnvnator_start)<=min(m.gridss_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and (abs(m.gridss_end-c.cnvnator_end)<=min(m.gridss_abs_svlen_20percent,c.cnvnator_abs_svlen_20percent)) and ((min(m.gridss_end,c.cnvnator_end)-min(m.gridss_start,c.cnvnator_start))>=max(m.gridss_abs_svlen_70percent,c.cnvnator_abs_svlen_70percent)))   )')

  # merge7_1 contains variables in gridss and cnvnator, regardless of whether they exist in manta, and some of their data columns

  merge7_2 = unique( merge7_1[,c("gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt", "cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size")] )
  # merge7_2 contains variables in gridss and cnvnator, regardless of whether they exist in manta, and only their key columns

  merge7_3 = sqldf('select * from merge7_2 t1 where not exists (select * from merge2 t2 where (t1.gridss_chrom=t2.gridss_chrom and t1.gridss_start=t2.gridss_start and t1.gridss_end=t2.gridss_end and t1.gridss_ref=t2.gridss_ref and t1.gridss_alt=t2.gridss_alt))')
  # merge7_3 contains variables in gridss and cnvnator that do not exist in manta, and only their key columns

  if (nrow(merge7_3) > 0) {

    we_have_gridss_cnvnator = "yes"

    merge7_3$manta_chrom = ''
    merge7_3$manta_start = 0
    merge7_3$manta_end = 0
    merge7_3$manta_ref = ''
    merge7_3$manta_alt = ''

    merge7_4 = merge7_3[,c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt", "gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt", "cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size")]
  }
}

# We now have the list of variables that appear in more than one data source (the data sources are manta, gridss, and cnvnator), with only their key fields. 
# Let's get all the data fields for all these variants.

if (we_have_cnvnator == "yes") {
  bind_all = merge5
} else {
  bind_all = merge5_2
}
if (we_have_manta_cnvnator == "yes") {
  bind2 = rbind( bind_all, merge6_4 )
  bind_all = bind2
}
if (we_have_gridss_cnvnator == "yes") {
  bind2 = rbind( bind_all, merge7_4 )
  bind_all = bind2
}

bind_all_manta = merge( x=bind_all, y=manta, by=c("manta_chrom", "manta_start", "manta_end", "manta_ref", "manta_alt"), all.x=TRUE, all.y=FALSE )
bind_all_manta[is.na(bind_all_manta)] = ""
# bind_all_manta contains variables in all 3 of manta and gridss and cnvnator or variables in only manta and gridss, and key columns of all 3 data sources, and manta data fields

bind_all_manta_gridss = merge( x=bind_all_manta, y=gridss, by=c("gridss_chrom", "gridss_start", "gridss_end", "gridss_ref", "gridss_alt"), all.x=TRUE, all.y=FALSE )
bind_all_manta_gridss[is.na(bind_all_manta_gridss)] = ""
# bind_all_manta contains variables in all 3 of manta and gridss and cnvnator or variables in only manta and gridss, and key columns of all 3 data sources, and manta and gridss data fields

if (we_have_cnvnator == "yes") {

  bind_all_manta_gridss_cnvnator = merge( x=bind_all_manta_gridss, y=cnvnator, by=c("cnvnator_chrom", "cnvnator_start", "cnvnator_end", "cnvnator_cnv_type", "cnvnator_cnv_size"), all.x=TRUE, all.y=FALSE )
  bind_all_manta_gridss_cnvnator[is.na(bind_all_manta_gridss_cnvnator)] = ""
  # bind_all_manta contains variables in all 3 of manta and gridss and cnvnator or variables in only manta and gridss, and key columns of all 3 data sources, and manta and gridss data fields

} else {
  bind_all_manta_gridss_cnvnator = bind_all_manta_gridss
}

bind_all_manta_gridss_cnvnator$manta_start = ifelse( (bind_all_manta_gridss_cnvnator$manta_start==0), '', bind_all_manta_gridss_cnvnator$manta_start )
bind_all_manta_gridss_cnvnator$manta_end = ifelse( (bind_all_manta_gridss_cnvnator$manta_end==0), '', bind_all_manta_gridss_cnvnator$manta_end )
bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent )
bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent )
bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent )

bind_all_manta_gridss_cnvnator$gridss_start = ifelse( (bind_all_manta_gridss_cnvnator$gridss_start==0), '', bind_all_manta_gridss_cnvnator$gridss_start )
bind_all_manta_gridss_cnvnator$gridss_end = ifelse( (bind_all_manta_gridss_cnvnator$gridss_end==0), '', bind_all_manta_gridss_cnvnator$gridss_end )
bind_all_manta_gridss_cnvnator$gridss_abs_svlen_15percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent )
bind_all_manta_gridss_cnvnator$gridss_abs_svlen_20percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent )
bind_all_manta_gridss_cnvnator$gridss_abs_svlen_70percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent )

if (we_have_cnvnator == "yes") {

  bind_all_manta_gridss_cnvnator$cnvnator_start = ifelse( (bind_all_manta_gridss_cnvnator$cnvnator_start==0), '', bind_all_manta_gridss_cnvnator$cnvnator_start )
  bind_all_manta_gridss_cnvnator$cnvnator_end = ifelse( (bind_all_manta_gridss_cnvnator$cnvnator_end==0), '', bind_all_manta_gridss_cnvnator$cnvnator_end )
  bind_all_manta_gridss_cnvnator$cnvnator_cnv_size = ifelse( (bind_all_manta_gridss_cnvnator$cnvnator_cnv_size==0), '', bind_all_manta_gridss_cnvnator$cnvnator_cnv_size )

  bind_all_manta_gridss_cnvnator$gridss_start = ifelse( (bind_all_manta_gridss_cnvnator$gridss_start==0), '', bind_all_manta_gridss_cnvnator$gridss_start )
  bind_all_manta_gridss_cnvnator$gridss_end = ifelse( (bind_all_manta_gridss_cnvnator$gridss_end==0), '', bind_all_manta_gridss_cnvnator$gridss_end )
  bind_all_manta_gridss_cnvnator$gridss_abs_svlen_15percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_15percent )
  bind_all_manta_gridss_cnvnator$gridss_abs_svlen_20percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_20percent )
  bind_all_manta_gridss_cnvnator$gridss_abs_svlen_70percent = ifelse( (bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent==999999999999), '', bind_all_manta_gridss_cnvnator$manta_abs_svlen_70percent )
}

bind_all = bind_all_manta_gridss_cnvnator

bind_all$sort_chrom = bind_all$manta_chrom
bind_all$sort_start = bind_all$manta_start
bind_all$sort_end = bind_all$manta_end
bind_all$sort_ref = bind_all$manta_ref
bind_all$sort_alt = bind_all$manta_alt
bind_all$sort_chrom = ifelse( (bind_all$sort_chrom==''), bind_all$gridss_chrom, bind_all$sort_chrom )
bind_all$sort_start = ifelse( (bind_all$sort_start==0), bind_all$gridss_start, bind_all$sort_start )
bind_all$sort_end = ifelse( (bind_all$sort_end==0), bind_all$gridss_end, bind_all$sort_end )
bind_all$sort_ref = ifelse( (bind_all$sort_ref==''), bind_all$gridss_ref, bind_all$sort_ref )
bind_all$sort_alt = ifelse( (bind_all$sort_alt==''), bind_all$gridss_alt, bind_all$sort_alt )

bind_all_sorted = bind_all[with(bind_all, order(sort_chrom, sort_start, sort_end, sort_ref, sort_alt)), ]
bind_all_sorted$sort_chrom = NULL
bind_all_sorted$sort_start = NULL
bind_all_sorted$sort_end = NULL
bind_all_sorted$sort_ref = NULL
bind_all_sorted$sort_alt = NULL

out_all = bind_all_sorted

out_all[is.na(out_all)] = ""
out_all[(out_all==".")] = ""

out_all$manta_abs_svlen_15percent = NULL
out_all$manta_abs_svlen_20percent = NULL
out_all$manta_abs_svlen_70percent = NULL
out_all$gridss_abs_svlen_15percent = NULL
out_all$gridss_abs_svlen_20percent = NULL
out_all$gridss_abs_svlen_70percent = NULL
if (we_have_cnvnator == "yes") {
  out_all$cnvnator_abs_svlen_15percent = NULL
  out_all$cnvnator_abs_svlen_20percent = NULL
  out_all$cnvnator_abs_svlen_70percent = NULL
}

write.table( out_all, file=outfile_all, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


