#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
indir_hrun=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt
in_ped_file="${batch}".pedigree_file.ped

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"
echo 'out_manifest:' $out_manifest

for infile in "${indir}"/*.CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.hot.compareSamples.tsv; do

  echo 'infile' $infile
  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  family_and_samples="${array[0]}"
  IFS='_' read -r -a array <<< "$family_and_samples"
  family_or_sample="${array[0]}" # this is either the sample id or it is the family id if there is one
  IFS='-' read -r -a array <<< "$family_or_sample"
  family_or_sample="${array[0]}" # this is either the sample id or it is the DCM family id if there is one
  echo 'family_or_sample' $family_or_sample 'in_ped_file' $in_ped_file
  samples_proband=`grep -P "^${family_or_sample}\t" "${in_ped_file}" | head -n 1`
  echo 'samples_proband' $samples_proband
  IFS=$'\t' read -r -a array <<< "$samples_proband"
  sample="${array[1]}"
  echo 'sample' $sample

  infile_tsv=$infile
  outfile_basename_prefix="${family_and_samples}".manta

  #samples_cohort=`grep -P "\t${family_or_sample}$|\t${sample}\t" "${in_manifest}"`
  samples_cohort=`grep -P "\t${family_or_sample}$|\t${sample}\t|\t${sample}$" "${in_manifest}"`
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"

  large_homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.large_homozygosity_runs_only.txt
  homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.txt

  echo -e "${family_and_samples}\t${infile_tsv}\t${outdir}\t${outfile_basename_prefix}\t${cohort}\t${large_homozygous_runs}\t${homozygous_runs}" >> "${out_manifest}"
  echo ''

done

echo 'out_manifest:' $out_manifest

./strucVars_extract_exon_hits_and_add_flags_SubmitJobs.sh $out_manifest

