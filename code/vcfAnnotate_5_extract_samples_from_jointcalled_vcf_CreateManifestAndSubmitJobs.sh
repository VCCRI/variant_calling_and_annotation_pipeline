#!/bin/bash
set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

infile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_annotate/"${batch}".gatk_vqsr_and_platypus_mnp.AllChrom.annovar_clinvarDATE.vep_reformatted.extras.vcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_annotate_samples

module load bcftools

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for sample in `bcftools query -l $infile`; do

  outfile_basename="${sample}".gatk_vqsr_and_platypus_mnp.AllChrom.annovar_clinvarDATE.vep_reformatted.extras.vcf.gz

  echo -e "${sample}\t${infile}\t${outdir}\t${outfile_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./extract_sample_from_jointcalled_vcf_SubmitJobs.sh $out_manifest

