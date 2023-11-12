#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile_done in "${indir}"/*.markDup.setTags.bqsr.bam.done; do

  infile="${infile_done%.done}"
  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./cnvnator_1_callCNVs_SubmitJobs.sh $out_manifest

