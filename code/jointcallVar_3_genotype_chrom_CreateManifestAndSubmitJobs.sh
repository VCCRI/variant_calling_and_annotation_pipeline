#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf_GenomicsDB # ./chr20
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf_GenomicsDB_genotype_vcf # ./chr20.vcf.gz

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
  infile="${indir}"/"${chrom}"

  echo -e "${chrom}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest
echo ''

./jointcallVar_3_genotype_chrom_SubmitJobs.sh $out_manifest

