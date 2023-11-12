#!/bin/bash

genes_without_bars=$1

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

reference=$exons_chrom_start_end_gene_strand

readarray -t array < $genes_without_bars

for gene in "${array[@]}"; do
  grep_string="\t${gene}\t"
  grep_string2=${grep_string//\./\\\.}
  grep_result=`grep -P "${grep_string2}" $reference`
  if [[ $grep_result == "" ]]; then
    echo 'Gene is not in reference:' $gene
  fi
done






