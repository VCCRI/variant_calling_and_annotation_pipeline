#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP
# SAMPLE_ONE.vqsr.decompose_blocksub_normalize.chrom1.MNP.vcf.gz
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

chr_prefix="chr"
if [[ $genome_version == "hg19" ]]; then
  chr_prefix=""
fi
#for infile in "${indir}"/*.vqsr.decompose_blocksub_normalize.chrom1.MNP.vcf.gz; do
for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp."${chr_prefix}"1.removeAltAsterisk.vcf; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  #infile_prefix="${indir}"/"${sample}".vqsr.decompose_blocksub_normalize.chrom
  #infile_suffix=".MNP.vcf.gz"
  #outfile="${outdir}"/"${sample}".AllChrom.MNP.vcf

  infile_prefix="${indir}"/"${sample}".gatk_vqsr_and_platypus_mnp.chr
  if [[ $genome_version == "hg19" ]]; then
    infile_prefix="${indir}"/"${sample}".gatk_vqsr_and_platypus_mnp.
  fi
  infile_suffix=.removeAltAsterisk.vcf
  outfile="${outdir}"/"${sample}".gatk_vqsr_and_platypus_mnp.AllChrom.removeAltAsterisk.vcf

  echo -e "${sample}\t${infile_prefix}\t${infile_suffix}\t${outdir}\t${outfile}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./concat_vcf_gz_chroms_SubmitJobs.sh $out_manifest

