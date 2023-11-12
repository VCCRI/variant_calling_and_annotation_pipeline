#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP # MY_BATCH.gatk_vqsr_and_platypus_mnp.chr22.vcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate # MY_BATCH.gatk_vqsr_and_platypus_mnp.chr10.annovar_clinvarDATE.hg38_multianno.vcf

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp.*.removeAltAsterisk.vcf; do

    infile_basename=$(basename $infile)
    IFS='.' read -r -a array <<< "$infile_basename"
    batch="${array[0]}"

    echo -e "${batch}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./vcfAnnotate_1_annotate_jointcalled_chrom_annovar_clinvar_SubmitJobs.sh $out_manifest

