#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate
tool=priority
param_file=/my/directory/variant_calling_and_annotation_pipeline/code/default_VPOT_filtered_for_SV.txt

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt
in_ped_file="${batch}".pedigree_file.ped

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"
echo 'out_manifest:' $out_manifest

for infile in "${indir}"/*_genes.tsv; do

  echo 'infile' $infile
  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  family_and_samples="${array[0]}"
  IFS='_' read -r -a array <<< "$family_and_samples"
  family_or_sample="${array[0]}" # this is either the sample id or it is the family id if there is one
  len=${#array[@]}
  family=''
  sample=''
  if [[ $len -eq 1 ]]; then
    is_singleton=1
    is_family=0
    sample=$family_or_sample # =$family_and_samples
    string_for_grep="\t${sample}\t"
  else
    is_family=1
    is_singleton=0
    family=$family_or_sample
    string_for_grep="^${family}\t"
  fi

  samples_cohort=`grep -P "${string_for_grep}" "${in_ped_file}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  infile_basename=$(basename $infile)
  outfile_prefix_basename="${infile_basename%.tsv}".vpot

  # Call vpot for singletons, and also for famiies.
  # Merged singletons report will get vpot-prioritised variants with vpot score > 0.
  # Merged families report will get all vpot-priroitised variants because we will do a manual inheritance modelling.

  echo -e "${family_and_samples}\t${tool}\t${param_file}\t${infile}\t${outdir}\t${outfile_prefix_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./call_vpot_sv_priority_SubmitJobs.sh $out_manifest

