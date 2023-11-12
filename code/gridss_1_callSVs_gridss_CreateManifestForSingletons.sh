#!/bin/bash
set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss
pedigree_file=/my/directory/variant_calling_and_annotation_pipeline/code/"${batch}"_ped_file.ped

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest.gridss.single_samples.txt

:>"${out_manifest}"

while read Family Subject Father Mother Sex Phenotype no_more; do
  if [[ $Family == $Subject ]]; then
    sample=$Subject
    infile="${indir}"/"${sample}".markDup.setTags.bqsr.bam
    echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"
  fi
done < $pedigree_file

echo ''
echo 'Run this next: ./gridss_1_callSVs_gridss_SubmitJobs.sh' $out_manifest
echo ''
