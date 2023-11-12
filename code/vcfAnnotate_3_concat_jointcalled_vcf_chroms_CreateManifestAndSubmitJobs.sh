#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_annotate

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

infiles="${indir}"/*.gatk_vqsr_and_platypus_mnp.chr1.annovar_clinvarDATE.vep_reformatted.extras.vcf
if [[ $genome_version == "hg19" ]]; then
  infiles="${indir}"/*.gatk_vqsr_and_platypus_mnp.1.annovar_clinvarDATE.vep_reformatted.extras.vcf
fi

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in $infiles; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  cohort="${array[0]}"

  infile_prefix_dirname=$(dirname $infile)
  infile_prefix="${infile_prefix_dirname}"/"${cohort}".gatk_vqsr_and_platypus_mnp.chr
  if [[ $genome_version == "hg19" ]]; then
    infile_prefix="${infile_prefix_dirname}"/"${cohort}".gatk_vqsr_and_platypus_mnp.
  fi
  infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.vcf
  outfile_basename="${cohort}".gatk_vqsr_and_platypus_mnp.AllChrom.annovar_clinvarDATE.vep_reformatted.extras.vcf
  echo -e "${cohort}\t${infile_prefix}\t${infile_suffix}\t${outdir}\t${outfile_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./concat_vcf_chroms_SubmitJobs.sh $out_manifest

