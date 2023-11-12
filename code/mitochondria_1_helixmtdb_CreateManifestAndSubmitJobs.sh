#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_MNP_annotate_samples
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/mitochondria

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

#for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp.AllChrom.annovar_clinvarDATE.vep_reformatted.extras.vcf.gz; do
for infile in "${indir}"/*.annovar_clinvarDATE_spliceai.vep.spliceogen.homrun.vcf; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}\t${outdir}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./mitochondria_1_helixmtdb_SubmitJobs.sh $out_manifest

