#!/bin/bash
set -euo pipefail

indir1=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_annotate_samples
indir2=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile1 in "${indir1}"/*.gatk_vqsr_and_platypus_mnp.AllChrom.annovar_clinvarDATE.vep_reformatted.extras.vcf.gz; do

  sample=$(basename $infile1)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"
  IFS='-' read -r -a array <<< "$sample"
  sample="${array[0]}"
  infile2="${indir2}"/output/"${sample}".gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz_out.txt

  echo -e "${sample}\t${infile1}\t${infile2}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./spliceogen_spliceai_2_annotate_vcf_with_spliceogen_SubmitJobs.sh $out_manifest

