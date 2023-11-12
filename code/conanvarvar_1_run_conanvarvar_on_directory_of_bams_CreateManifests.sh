#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/conanvarvar

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

i=0
num_samples_per_conanvarvar_batch=10
for infile in "${indir}"/*.bam; do

  i=$(( i + 1 ))
  div_i=$(( $i / $num_samples_per_conanvarvar_batch ))
  mod_i=$(( $i % $num_samples_per_conanvarvar_batch ))
  if [[ "$mod_i" != 0 ]]; then
    div_i=$(( $div_i + 1 ))
  fi
  if [[ "$mod_i" == 1 ]]; then
    newdir="${outdir}/group_${div_i}"
    mkdir -p $newdir
    new_manifest=./manifests/manifest."${curr_pgm_base}".group_"${div_i}".txt
    :>$new_manifest
    echo -e "${newdir}\t${newdir}" >> $new_manifest
    #echo 'outfile manifest' $div_i ':' $new_manifest
    echo '# The directories and file links have been set up. Please run the following command:'
    echo './conanvarvar_1_run_conanvarvar_on_directory_of_bams_SubmitJobs.sh' $new_manifest
    echo ''
  fi
  infile_basename=$(basename $infile)
  new_link="${newdir}"/"${infile_basename}"
  infile_bai="${infile%.bam}".bai
  infile_bai_basename=$(basename $infile_bai)
  infile_bai_basename_modified=${infile_bai_basename/.bam/}
  new_link_bai="${newdir}"/"${infile_bai_basename_modified}"
  if [[ -h "$new_link" ]]; then
    :
  else
    ln -s $infile $newdir
  fi
  if [[ -h "$new_link_bai" ]]; then
    :
  else
    ln -s $infile_bai $new_link_bai
  fi
done


