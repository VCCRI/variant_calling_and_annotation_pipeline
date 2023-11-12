#!/bin/bash

# This script extracts only the first reference record of regions of a gene. Usuallyl there is only one per gene.

module load bedtools

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

# Input:
genes_without_bars=$1
reference=$ucsc_refseq_genes

# Output:
genes_basename=$(basename $genes_without_bars)
genes_dirname=$(dirname $genes_without_bars)
genes_basename="${genes_basename%.txt}"
genes_basename="${genes_basename%.tsv}"
genes_basename="${genes_basename%.bed}"
genes_basename="${genes_basename%_withoutBars}"
genes_with_bars="${genes_dirname}"/"${genes_basename}"_withBars.txt
regions_tab_delimited="${genes_dirname}"/"${genes_basename}"_"${genome_version}"_RefSeq_regions_tab_delimited.txt
regions_colon_dash_format="${genes_dirname}"/"${genes_basename}"_"${genome_version}"_RefSeq_regions_colon_dash_format.txt

temp_regions_tab_delimited='temp_regions_tab_delimited.bed'

awk 'BEGIN {FS=" ";OFS=""} {print "|" $1 "|"}' $genes_without_bars > $genes_with_bars

readarray -t array < $genes_without_bars

:>$temp_regions_tab_delimited
for gene in "${array[@]}"; do
  grep_string="\t${gene}\t"
  grep_result=`grep -P "${grep_string}" $reference`
  if [[ $grep_result != "" ]]; then
    IFS=$'\t' read -r -a array2 <<< "$grep_result"
    chrom="${array2[0]}"
    start="${array2[1]}"
    end="${array2[2]}"
    out_string1="${chrom}:${start}-${end}"
    out_string2="${chrom}\t${start}\t${end}"
    echo -e "${out_string2}" >> $temp_regions_tab_delimited
  fi
done

bedtools sort -i $temp_regions_tab_delimited | bedtools merge -i - > $regions_tab_delimited

awk 'BEGIN {FS="\t";OFS=""} {print $1 ":" $2 "-" $3}' $regions_tab_delimited > $regions_colon_dash_format


