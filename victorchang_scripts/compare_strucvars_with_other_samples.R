
options(width=200)
options(stringsAsFactors = FALSE)
library(tidyr)
library(dplyr)
library(data.table)
library(sqldf)

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( 'A020029517111', 'A020029517111.manta.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.bndInsideExon.tsv', 'A020029517111.manta.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.bndInsideExon.compareSamples.tsv', 'database.manta.allVars.txt', 'sample1|sample2|sample3' )

sample = as.character(args[1])
infile = as.character(args[2])
outfile = as.character(args[3])
database_file = as.character(args[4])
include_low_qual_flag = as.character(args[5])

exclude_samples = list()
if (length(args) > 5) {
  exclude_samples_list = as.character(args[6])
  exclude_samples = unlist(strsplit(exclude_samples_list, '\\|'))
  print('exclude_samples_list')
  print(exclude_samples_list)
}

print('exclude_samples')
print(exclude_samples)

indata = read.table( infile, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
database = read.table( database_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(database) = c("sample","filename")

# sqldf crashes when got same_name, different_capitals colnames
if ( ("end" %in% colnames(indata)) & ("END" %in% colnames(indata)) ) {
  names(indata)[names(indata) == 'END'] = 'INFO_END'
}
if ( ("End" %in% colnames(indata)) & ("END" %in% colnames(indata)) ) {
  names(indata)[names(indata) == 'END'] = 'INFO_END'
}
# sqldf crashes when got same_name, different_capitals colnames
if ( ("ref" %in% colnames(indata)) & ("REF" %in% colnames(indata)) ) {
  names(indata)[names(indata) == 'REF'] = 'INFO_REF'
}
if ( ("Ref" %in% colnames(indata)) & ("REF" %in% colnames(indata)) ) {
  names(indata)[names(indata) == 'REF'] = 'INFO_REF'
}

fill_svlen_column <- function( indata ) {

  colnames(indata)[1] = "chrom"
  colnames(indata)[2] = "start"
  colnames(indata)[3] = "end"
  names(indata)[names(indata) == 'SVLEN'] = 'svlen'
  indata$chrom = as.character(indata$chrom)
  indata$start = as.numeric(indata$start)
  indata$end = as.numeric(indata$end)

  if ("svlen" %in% colnames(indata)) {
    do_nothing = 1
  } else {
    if ("cnvlen" %in% colnames(indata)) {
      names(indata)[names(indata) == "cnvlen"] = "svlen"
    } else {
      indata$svlen = indata$end - indata$start + 1
    }
  }
  indata$svlen = as.numeric(indata$svlen)
  indata$svlen = ifelse( is.na(indata$svlen), (indata$end - indata$start + 1), indata$svlen )
  indata$svlen = abs(indata$svlen)
  indata
}

if (nrow(indata) > 0) {

  indata = fill_svlen_column( indata )

  indata$svlen_5percent = round( indata$svlen * 0.05, digits=0 )

  save_colname_4 = colnames(indata)[4]
  save_colname_5 = colnames(indata)[5]
  colnames(indata)[4] = "ref"
  colnames(indata)[5] = "alt"
  indata$ref = as.character(indata$ref)
  indata$alt = as.character(indata$alt)

  indata$newCol_otherSamples_exactMatch = ""
  indata$newCol_otherSamples_90percentMatch = ""

  for ( i in 1:nrow(database) ) {

    this_row = database[i,]
    otherSample_sample = as.character(this_row$sample)
    otherSample_filename = as.character(this_row$filename)

    otherSample_is_an_excluded_sample = FALSE
    if (otherSample_sample == sample) {
      otherSample_is_an_excluded_sample = TRUE
    }
    if (length(exclude_samples) > 0) {
      for ( j in 1:length(exclude_samples) ) {
        this_excluded_sample = exclude_samples[j]
        if (otherSample_sample == this_excluded_sample) {
          otherSample_is_an_excluded_sample = TRUE
        }
      }
    }

    if (!otherSample_is_an_excluded_sample) {

      print('otherSample_filename')
      print(otherSample_filename)
      file_info = file.info(otherSample_filename)
      print('file_info')
      print(file_info)
      is_file_empty = file_info$size

      if (is_file_empty != 0) {

        #otherSample_sample = "19W001098"
        #otherSample_filename = "19W001098.gridss_1_4_1.CDSexons.EHRFr99.gnomadSV.DGV.SegDup.dbscSNV.entireGene.collapsed.CHD_tier1_106_genes_plus_extended_505_genes_plus_CFAP54.tsv"
        otherSample = read.table( otherSample_filename, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )

        if (nrow(otherSample) > 0) {

          otherSample = fill_svlen_column( otherSample )

          colnames(otherSample)[1] = "chrom2"
          colnames(otherSample)[2] = "start2"
          colnames(otherSample)[3] = "end2"
          colnames(otherSample)[4] = "ref2"
          colnames(otherSample)[5] = "alt2"
          otherSample$ref2 = as.character(otherSample$ref2)
          otherSample$alt2 = as.character(otherSample$alt2)
          otherSample$sample = otherSample_sample

          otherSample1 = otherSample[,c("sample","chrom2","start2","end2","ref2","alt2","svlen")]
          otherSample1$svlen_5percent_2 = round( otherSample1$svlen * 0.05, digits=0 )
          otherSample1 = unique(otherSample1)

          # sqldf crashes when got same_name, different_capitals colnames
          if ( ("end" %in% colnames(otherSample1)) & ("END" %in% colnames(otherSample1)) ) {
            names(otherSample1)[names(otherSample1) == 'END'] = 'END_COLUMN'
          }

          merge1 = sqldf("select * from indata d, otherSample1 s where d.chrom=s.chrom2 and d.start=s.start2 and d.end=s.end2")
          merge1 = unique(merge1)

          merge2 = merge1[,c("chrom","start","end","ref","alt","sample")]
          merge2 = unique(merge2)

          merge3 = merge( x=indata, y=merge2, by=c("chrom","start","end","ref","alt"), all.x=TRUE, all.y=FALSE )
          merge3[is.na(merge3)] = ""

          merge3$newCol_otherSamples_exactMatch = ifelse( (merge3$newCol_otherSamples_exactMatch==""), merge3$sample, paste(merge3$newCol_otherSamples_exactMatch,merge3$sample,sep=", ") )
          merge3$newCol_otherSamples_exactMatch = gsub( ", $", "", merge3$newCol_otherSamples_exactMatch )

          merge3$sample = NULL
          merge3 = unique(merge3)

          indata = merge3

          rm(merge1,merge2,merge3)

          merge1 = sqldf("select * from indata d, otherSample1 s where (d.chrom=s.chrom2) and (abs(d.start-s.start2)<=min(d.svlen_5percent,s.svlen_5percent_2)) and (abs(d.end-s.end2)<=min(d.svlen_5percent,s.svlen_5percent_2)) and ((min(d.end,s.end2)-min(d.start,s.start2))>=max(d.svlen_5percent,s.svlen_5percent_2))")
          merge1 = unique(merge1)

          merge2 = merge1[,c("chrom","start","end","ref","alt","sample")]
          merge2 = unique(merge2)

          merge3 = merge( x=indata, y=merge2, by=c("chrom","start","end","ref","alt"), all.x=TRUE, all.y=FALSE )
          merge3[is.na(merge3)] = ""

          merge3$newCol_otherSamples_90percentMatch = ifelse( (merge3$newCol_otherSamples_90percentMatch==""), merge3$sample, paste(merge3$newCol_otherSamples_90percentMatch,merge3$sample,sep=", ") )
          merge3$newCol_otherSamples_90percentMatch = gsub( ", $", "", merge3$newCol_otherSamples_90percentMatch )

          merge3$sample = NULL
          merge3 = unique(merge3)

          indata = merge3

          rm(merge1,merge2,merge3)
        }
      }
    }
  }

  indata$svlen_5percent = NULL

  colnames(indata)[4] = save_colname_4
  colnames(indata)[5] = save_colname_5
}

if ((include_low_qual_flag == '') | (include_low_qual_flag == '.')) {
  names(indata)[names(indata) == 'newCol_otherSamples_exactMatch'] = 'otherSamples_exactMatch'
  names(indata)[names(indata) == 'newCol_otherSamples_90percentMatch'] = 'otherSamples_90percentMatch'
} else {
  names(indata)[names(indata) == 'newCol_otherSamples_exactMatch'] = 'includeLowQual_otherSamples_exactMatch'
  names(indata)[names(indata) == 'newCol_otherSamples_90percentMatch'] = 'includeLowQual_otherSamples_90percentMatch'
}
write.table( indata, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


