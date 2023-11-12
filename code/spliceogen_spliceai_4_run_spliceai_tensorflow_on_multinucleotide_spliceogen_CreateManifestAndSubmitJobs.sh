#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai
indir2=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP/*.vcf # just to get a good vcf header with all the chromosomes
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

infile2=`ls -1 $indir2 | head -n 1`

:>"${out_manifest}"

for infile in "${indir}"/*.needs_spliceai_tensorflow_scores.tsv; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}\t${infile2}\t${outdir}" >> "${out_manifest}"
done

echo 'out_manifest:' $out_manifest

./spliceogen_spliceai_4_run_spliceai_tensorflow_on_multinucleotide_spliceogen_SubmitJobs.sh $out_manifest

