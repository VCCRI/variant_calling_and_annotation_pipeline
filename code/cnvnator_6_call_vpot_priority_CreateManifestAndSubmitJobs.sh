#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate
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

# SAMPLE.cnvnator.CHD_955_genes.tsv
# SAMPLE.cnvnator.DCM_909_genes.tsv
for infile in "${indir}"/*_genes.tsv; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"
  sample_for_grep=$sample

  sample_for_grep="\t${sample_for_grep}\t|\t${sample_for_grep}$"
  samples_cohort=`grep -P "${sample_for_grep}" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  infile_basename=$(basename $infile)
  outfile_prefix_basename="${infile_basename%.tsv}".vpot

  echo -e "${sample}\t${tool}\t${param_file}\t${infile}\t${outdir}\t${outfile_prefix_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./call_vpot_sv_priority_SubmitJobs.sh $out_manifest

