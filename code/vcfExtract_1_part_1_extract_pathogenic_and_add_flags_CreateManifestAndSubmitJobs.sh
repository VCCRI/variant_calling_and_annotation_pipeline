#!/bin/bash
set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples
infile_prefix_substring=.gatk_vqsr_and_platypus_mnp.chr
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.homrun.vcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples_extracts
if [[ $genome_version == "hg19" ]]; then
  infile_prefix_substring=.gatk_vqsr_and_platypus_mnp.
fi

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*"${infile_prefix_substring}"1"${infile_suffix}"; do

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

./vcfExtract_1_part_1_extract_pathogenic_and_add_flags_SubmitJobs.sh $out_manifest

