#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

manifest="${batch}".one_line_per_family.txt
indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples_extracts
infile_basename_suffix=.gatk_platypus_mnp.AllChrom.annovar_vep.stopLossGain.vpot_final_output_file.vpot_gt_0.txt
outfile_basename_suffix=_"${batch}"__gatk_platypus_mnp.AllChrom.annovar_vep.stopLossGain.vpot_final_output_file.vpot_gt_0.txt

tmphdr='temp_hdr.txt'
tmphdr2='temp_hdr2.txt'

cohorts=(CHD DCM)

for cohort in "${cohorts[@]}"; do
  echo 'Processing:' $cohort
  outfile="${indir}"/"${cohort}${outfile_basename_suffix}"
  string_of_files=''
  while IFS= read -r inline; do
    IFS=$'\t' read -r -a array <<< "$inline"
    this_flagship="${array[1]}"
    this_sample="${array[3]}"
    this_cohort="${array[0]}"
    if [[ "$this_cohort" == "$cohort" ]]; then
      infile="${indir}"/"${this_sample}${infile_basename_suffix}"
      new_string_of_files="${string_of_files} ${infile}"
      string_of_files=$new_string_of_files
      one_file=$infile
    fi
  done < "$manifest"
  head -n 1 $one_file > $tmphdr
  last_col=$(sed -e 's/\t/\n/g' "$tmphdr" | tail -n 1)
  sed_string='s/'"${last_col}"'$/sample_field/'
  sed $sed_string $tmphdr > $tmphdr2
  cat $string_of_files | grep -P -v '^Rank|^chrom' | sort -k1nr -k2nr | cat $tmphdr2 - > $outfile || true
  echo ''
  echo 'outfile:' $outfile
  echo ''

done
