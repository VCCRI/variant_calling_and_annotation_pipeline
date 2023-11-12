#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate

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
echo 'out_manifest:' $out_manifest

for infile in "${indir}"/*.manta_diploidSV_withInversions.rmvCols.tsv; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"
  IFS='_' read -r -a array <<< "$sample"
  family="${array[0]}"
  infile_basename=$(basename $infile)
  outfile_basename="${infile_basename%.txt}"
  outfile_basename="${outfile_basename%.tsv}"
  outfile_basename="${outfile_basename}".annotate.tsv

  echo 'family:' $family
  family_for_grep="\t${family}\t|\t${family}$"
  samples_cohort=`grep -P "${family_for_grep}" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  echo -e "${sample}\t${infile}\t${outdir}\t${cohort}" >> "${out_manifest}"
done

echo 'out_manifest:' $out_manifest

./manta_2_annotate_multiple_bedtools_hits_SubmitJobs.sh $out_manifest

