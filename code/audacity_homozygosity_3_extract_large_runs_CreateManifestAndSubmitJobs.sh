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

for infile in "${indir}"/*.audacity_homozygosity_output.txt; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./audacity_homozygosity_3_extract_large_runs_SubmitJobs.sh $out_manifest

