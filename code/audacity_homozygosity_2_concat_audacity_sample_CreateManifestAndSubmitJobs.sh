#!/bin/bash
set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*_chr1/*/*_DIDOH3M2Regions.txt; do

  infile_bed="${infile}".bed
  echo ' ln -s' $infile $infile_bed
  ln -s $infile $infile_bed
  echo ''

  infile_basename=$(basename $infile_bed)
  IFS='_' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"
  IFS='-' read -r -a array <<< "$sample"
  sample="${array[0]}"

  sample_for_grep="\t${sample}\t|\t${sample}$"
  echo 'samples_cohort=grep -P' "${sample_for_grep}" "${in_manifest}"
  echo ''
  samples_cohort=`grep -P "${sample_for_grep}" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  infile_prefix="${indir}/${sample}_chr"
  infile_suffix="/${sample}/${sample}_DIDOH3M2Regions.txt"

  echo -e "${sample}\t${infile_prefix}\t${infile_suffix}\t${outdir}\t${cohort}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./audacity_homozygosity_2_concat_audacity_sample_SubmitJobs.sh $out_manifest

