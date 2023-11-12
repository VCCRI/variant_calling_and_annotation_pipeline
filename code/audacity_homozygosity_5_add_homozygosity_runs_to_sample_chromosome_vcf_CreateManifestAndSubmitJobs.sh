#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples
indir_hrun=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.gatk_vqsr_and_platypus_mnp.chr*.annovar_clinvarDATE.vep_reformatted.extras.vcf; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"

  flag1=BigHomozygRun
  large_homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.large_homozygosity_runs_only.txt
  flag2=HRun
  homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.txt

  sample_job_differentiating_id="${array[2]}" # chromosome

  echo -e "${sample}\t${infile}\t${outdir}\t${flag1}\t${large_homozygous_runs}\t${flag2}\t${homozygous_runs}\t${sample_job_differentiating_id}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./audacity_homozygosity_add_homozygosity_runs_to_sample_vcf_EntryPoint.sh $out_manifest

