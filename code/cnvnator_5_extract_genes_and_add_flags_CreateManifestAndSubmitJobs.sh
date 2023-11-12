#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate
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

for infile in "${indir}"/*.CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.hot.compareSamples.tsv; do

  echo 'infile' $infile
  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"
  echo 'sample' $sample 'in_manifest' $in_manifest

  infile_tsv=$infile
  outfile_basename_prefix="${sample}".cnvnator

  samples_cohort=`grep -P "\t${sample}\t|\t${sample}$" "${in_manifest}"`
  echo 'samples_cohort' $samples_cohort
  IFS=$'\t' read -r -a array <<< "$samples_cohort"
  cohort="${array[0]}"
  echo 'cohort' $cohort

  large_homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.large_homozygosity_runs_only.txt
  homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.txt

  echo -e "${sample}\t${infile_tsv}\t${outdir}\t${outfile_basename_prefix}\t${cohort}\t${large_homozygous_runs}\t${homozygous_runs}" >> "${out_manifest}"
  echo ''

done

echo 'out_manifest:' $out_manifest

./strucVars_extract_genes_and_add_flags_SubmitJobs.sh $out_manifest

