#!/bin/bash
set -euo pipefail

#indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP
# SAMPLE_ONE.vqsr.decompose_blocksub_normalize.chrom1.MNP.vcf.gz
indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_samples
# SAMPLE_ONE.gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen

# or use
# /my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP/MY_BATCH.gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

#for infile in "${indir}"/*.AllChrom.MNP.vcf.gz; do
for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf.gz; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./spliceogen_1_launch_multiple_samples_to_pbs_SubmitJobs.sh $out_manifest

