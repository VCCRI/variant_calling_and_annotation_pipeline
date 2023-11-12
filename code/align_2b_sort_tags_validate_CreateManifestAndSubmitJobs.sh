#!/bin/bash
set -euo pipefail

#NOTE - markDuplicates step already completed
# with output file <sample>.markDup.bam
indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/fastq_bam_markDup_setTags
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/fastq_bam_markDup_setTags

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

#NOTE - only runs on already processed samples
for infile in "${indir}"/*.markDup.bam; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"
  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./align_2b_sort_tags_validate_SubmitJobs.sh $out_manifest

