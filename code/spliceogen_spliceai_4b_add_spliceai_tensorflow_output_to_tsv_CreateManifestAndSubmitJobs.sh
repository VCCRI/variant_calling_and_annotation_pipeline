#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.annovar_clinvarDATE_spliceai.vep.spliceogen.gnomad_filtered_splicing_variants.*.needs_spliceai_tensorflow_scores.tensorflow_results.tsv; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  infile2="${infile%.tensorflow_results.tsv}"
  infile2="${infile2%.needs_spliceai_tensorflow_scores}"
  infile2="${infile2}".needs_spliceai_tensorflow_scores.tsv

  echo -e "${sample}\t${infile}\t${infile2}\t${outdir}" >> "${out_manifest}"
done

echo 'out_manifest:' $out_manifest

./spliceogen_spliceai_4b_add_spliceai_tensorflow_output_to_tsv_SubmitJobs.sh $out_manifest

