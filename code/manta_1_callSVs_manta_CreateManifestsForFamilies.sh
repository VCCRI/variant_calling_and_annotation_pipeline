#!/bin/bash
set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta
outdir2=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
pedigree_file=/my/directory/variant_calling_and_annotation_pipeline/code/"${batch}"_ped_file.ped

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"

echo ''
echo 'Run the following next:'
prev_family=""
while read Family Subject Father Mother Sex Phenotype no_more; do
  if [[ $Family != "Family" ]]; then
    if [[ $Family != $Subject ]]; then
      if [[ $Family != $prev_family ]]; then
        out_manifest=./manifests/test_manifest.manta."${Family}".txt
        :>$out_manifest
        echo './manta_1_jointcallSVs_manta_SubmitJobs.sh' $out_manifest
        prev_family=$Family
      fi
      sample=$Subject
      infile="${indir}"/"${sample}".markDup.setTags.bqsr.bam
      echo -e "${sample}\t${Family}\t${infile}\t${outdir2}\t${outdir}" >> "${out_manifest}"
    fi
  fi
done < $pedigree_file
echo ''
