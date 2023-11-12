#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples_extracts
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples_extracts

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

ARG1=${1:-none}
if [ $ARG1 == "name" ]; then
  echo 'out_manifest (not created):' $out_manifest
  exit
fi

:>"${out_manifest}"

for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp.AllChrom.annovar_clinvarDATE.vep_reformatted.extras.homrun.stoploss_stopgain.tsv; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  sample_for_grep="\t${sample}\t|\t${sample}$"
  samples_cohort=`grep -P "${sample_for_grep}" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  tool='priority'

  param_file='/my/directory/variant_calling_and_annotation_pipeline/code/default_VPOT_MAF_0p01_PPF_for_CHD_INFOnames.txt'
  if [ "$cohort" == "DCM" ]; then
    param_file=/my/directory/variant_calling_and_annotation_pipeline/code/default_VPOT_MAF_0p001_PPF_for_DCM_INFOnames.txt
  fi

  out_prefix_basename="${sample}.gatk_platypus_mnp.AllChrom.annovar_vep.stopLossGain.vpot"

  echo -e "${sample}\t${tool}\t${param_file}\t${infile}\t${outdir}\t${out_prefix_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./call_vpot_vcf_tsv_single_samples_for_stopgain_SubmitJobs.sh $out_manifest

