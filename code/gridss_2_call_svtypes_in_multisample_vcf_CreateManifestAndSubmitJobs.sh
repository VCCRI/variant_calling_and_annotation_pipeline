#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.gridss.vcf.gz; do

  infile_basename=$(basename $infile)
  first3=$(echo $infile_basename | cut -c1-3)
  IFS='.' read -r -a array <<< "$infile_basename"
  family_and_samples="${array[0]}"
  IFS='_' read -r -a array <<< "$family_and_samples"
  array_length=${#array[@]}
  if [[ $array_length -gt 1 ]]; then
    family="${array[0]}"
    sample="${array[1]}"
    infile_basename=$(basename $infile)
    outfile_basename="${infile_basename%.sh}"
    outfile_basename="${outfile_basename%.gz}"
    outfile_basename="${outfile_basename%.vcf}"
    outfile_basename="${outfile_basename}".filter_svtypes.vcf
    echo -e "${family_and_samples}\t${infile}\t${outdir}\t${outfile_basename}" >> "${out_manifest}"
  fi

done

echo 'out_manifest:' $out_manifest

./gridss_2_call_svtypes_in_multisample_vcf_SubmitJobs.sh $out_manifest

