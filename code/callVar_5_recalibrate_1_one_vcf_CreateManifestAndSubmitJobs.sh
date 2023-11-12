#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf_GenomicsDB_genotype_vcf # MY_BATCH__allChrom.vcf.gz
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_vqsr # MY_BATCH__allChrom.vqsr_intermediate_file.*

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*_allChrom.vcf.gz; do

  infile_basename=$(basename $infile)
  new_string="${infile_basename/__/.}"
  IFS='.' read -r -a array <<< "$new_string"
  sample="${array[0]}"
  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./callVar_5_recalibrate_1_one_vcf_SubmitJobs.sh $out_manifest

