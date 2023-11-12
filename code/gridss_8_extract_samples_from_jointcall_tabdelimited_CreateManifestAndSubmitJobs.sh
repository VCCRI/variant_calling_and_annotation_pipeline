#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate
# FAMILY_FOUR_SAMPLE_ONE_SAMPLE_TWO_SAMPLE_THREE.gridss.CHD_all2022jun_955_genes.vpot_final_output_file.txt
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate
indir_hrun=/my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt
# CHD	A0133002	.	SAMPLE_ONE
# CHD	A1933066	.	SAMPLE_ONE
# CHD	A1433012	FAMILY_FOUR	SAMPLE_ONE	SAMPLE_TWO	SAMPLE_THREE
# CHD	A1133010	FAMILY_FIVE	SAMPLE_ONE	SAMPLE_TWO	SAMPLE_THREE

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

while IFS= read -r inline; do

  IFS=$'\t' read -r -a array <<< "$inline"
  cohort="${array[0]}"
  family="${array[2]}"
  if [[ $family != '.' ]]; then

    arr_len="${#array[@]}"
    arr_last=$(( $arr_len - 1 ))
    i=3
    while [[ $i -lt $arr_len ]] ; do

      sample="${array[$i]}"
      (( i += 1 ))

      #infile=$(ls -1 "${indir}"/"${family}"_*.gridss.*_genes.vpot_final_output_file.txt | head -n 1)
      infile=$(ls -1 "${indir}"/"${family}"*.gridss.*_genes.vpot_final_output_file.txt | head -n 1)
      infile_basename=$(basename $infile)
      IFS='.' read -r -a array2 <<< "$infile_basename"
      family_and_samples="${array2[0]}"
      outfile_basename="${infile_basename/$family_and_samples/$sample}"
      outfile="${outdir}"/"${outfile_basename}"

      large_homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.large_homozygosity_runs_only.txt
      homozygous_runs="${indir_hrun}"/"${sample}".audacity_homozygosity_output.txt

      echo -e "${sample}\t${infile}\t${outdir}\t${outfile}\t${large_homozygous_runs}\t${homozygous_runs}" >> "${out_manifest}"

    done
  fi
done < "$in_manifest"

echo 'out_manifest:' $out_manifest

./strucVars_extract_sample_from_jointcall_tabdelimited_SubmitJobs.sh $out_manifest

