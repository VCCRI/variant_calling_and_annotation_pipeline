#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf_GenomicsDB # chr4

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

:>"${out_manifest}"

mito="M"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
fi

for chrom_id in {1..22} "X" "Y" "$mito" ; do

  chrom="chr"$chrom_id
  if [[ $genome_version == "hg19" ]]; then
    chrom=$chrom_id
  fi
  infile_list=''

  for infile in "${indir}"/*.chrom"${chrom_id}".g.vcf.gz; do

    if [[ $infile_list == '' ]]; then
      infile_list=$infile
    else
      new_infile_list=$infile_list':'$infile
      infile_list=$new_infile_list
    fi

  done

  echo -e "${chrom}\t${outdir}\t${infile_list}" >> "${out_manifest}"

done

echo ''
echo 'out_manifest:' $out_manifest
echo ''
echo 'The GenomicsDB will be by chromosome, with multiple samples in each chromosome.'
echo ''

./jointcallVar_2_GenomicsDBImport_chrom_SubmitJobs.sh $out_manifest

