#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf_GenomicsDB_genotype_vcf # ./chr1.vcf.gz ./chr2.vcf.gz ./chr3.vcf.gz
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gvcf_GenomicsDB_genotype_vcf # ./MY_BATCH__allChrom.vcf.gz

currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

infile_prefix="${indir}"/chr
if [[ $genome_version == "hg19" ]]; then
  infile_prefix="${indir}"/
fi
infile_suffix=".vcf.gz"
outfile="${outdir}"/"${batch}"__allChrom.vcf.gz

echo -e "${batch}\t${infile_prefix}\t${infile_suffix}\t${outfile}" >> "${out_manifest}"

echo ''
echo 'out_manifest:' $out_manifest

./callVar_4_combine_chrom_SubmitJobs.sh $out_manifest

