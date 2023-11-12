#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai
indir_hrun=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.annovar_clinvarDATE_spliceai.vep.spliceogen.tsv; do # catch spliceogen and spliceai

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"
  sample_for_grep="\t${sample}\t|\t${sample}$"

  samples_cohort=`grep -P "${sample_for_grep}" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  large_homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.large_homozygosity_runs_only.txt
  homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.txt

  echo -e "${sample}\t${infile}\t${outdir}\t${cohort}\t${large_homozygous_runs}\t${homozygous_runs}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./spliceogen_spliceai_3_extract_genes_and_splice_variants_of_interest_SubmitJobs.sh $out_manifest

