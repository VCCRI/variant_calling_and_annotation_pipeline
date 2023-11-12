#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_vqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_vqsr_normalize

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*_allChrom.vqsr.vcf.gz; do

  infile_basename=$(basename $infile)
  new_string="${infile_basename/__/.}"
  IFS='.' read -r -a array <<< "$new_string"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./vcfAdjust_1_decompose_SubmitJobs.sh $out_manifest

