#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp.*.removeAltAsterisk.annovar_clinvarDATE."${genome_version}"_multianno.rmvSomeAnnovar.vep_reformatted.vcf; do

    infile_basename=$(basename $infile)
    IFS='.' read -r -a array <<< "$infile_basename"
    batch="${array[0]}"

    echo -e "${batch}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./vcfAnnotate_2b_annotate_jointcalled_chrom_extras_SubmitJobs.sh $out_manifest

