#!/bin/bash
set -euo pipefail

infile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP/MY_BATCH.gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_samples
# SAMPLE_ONE.gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz

module load bcftools

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for sample in `bcftools query -l $infile`; do

  outfile_basename="${sample}".gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz

  echo -e "${sample}\t${infile}\t${outdir}\t${outfile_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./extract_sample_from_jointcalled_vcf_SubmitJobs.sh $out_manifest

