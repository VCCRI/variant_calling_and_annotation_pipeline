#!/bin/bash

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

echo "sample" "total_hrun" "percent_hrun" "total_big_hrun" "percent_big_hrun"

for infile in "${indir}"/*.audacity_homozygosity_output.txt; do

  infile2="${infile%.txt}".large_homozygosity_runs_only.txt

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  result=$(awk -v sample="$sample" 'BEGIN {
    FS="\t"
    OFS=","
    sum_hrun = 0
  }
  {
  if (NR>1) {
    hrun_len = $3 - $2 + 1
    sum_hrun = sum_hrun + hrun_len
  }
  }
  END {
    percent=sum_hrun/30000000
    print sum_hrun, percent
  }' $infile)

  IFS=',' read -r -a array <<< "$result"
  hrun="${array[0]}"
  hrun_pc="${array[1]}"

  result=$(awk -v sample="$sample" 'BEGIN {
    FS="\t"
    OFS=","
    sum_hrun = 0
  }
  {
  if (NR>1) {
    hrun_len = $3 - $2 + 1
    sum_hrun = sum_hrun + hrun_len
  }
  }
  END {
    percent=sum_hrun/30000000
    print sum_hrun, percent
  }' $infile2)

  IFS=',' read -r -a array <<< "$result"
  big_hrun="${array[0]}"
  big_hrun_pc="${array[1]}"

  #declare -A array2
  #array2[908]=plot
  #array2[2035]=plot

  echo $sample $hrun $hrun_pc $big_hrun $big_hrun_pc # "${array2[$sample]}"

done




