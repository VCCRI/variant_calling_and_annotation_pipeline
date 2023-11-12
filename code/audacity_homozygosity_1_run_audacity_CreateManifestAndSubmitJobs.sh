#!/bin/bash
set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples
infile_prefix_substring=.gatk_vqsr_and_platypus_mnp.chr
if [[ $genome_version == "hg19" ]]; then
  infile_prefix_substring=.gatk_vqsr_and_platypus_mnp.
fi
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.vcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

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

  echo 'infile' $infile
  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  sample_for_grep="\t${sample}\t|\t${sample}$"
  samples_cohort=`grep -P "${sample_for_grep}" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  infile_prefix="${indir}/${sample}${infile_prefix_substring}"

  echo -e "${sample}\t${infile_prefix}\t${infile_suffix}\t${outdir}\t${cohort}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./audacity_homozygosity_1_run_audacity_SubmitJobs.sh $out_manifest

