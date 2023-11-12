#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.cnvnator.txt; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"
  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./cnvnator_2_filter_SubmitJobs.sh $out_manifest

