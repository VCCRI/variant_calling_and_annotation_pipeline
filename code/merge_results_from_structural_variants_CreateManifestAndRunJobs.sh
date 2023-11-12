#!/bin/bash
set -euo pipefail

indir_gridss=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate
# SAMPLE_THREE.gridss.CHD_955_genes.vpot_final_output_file.txt
indir_manta=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
indir_cnvnator=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate

outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/merge_results

infix_CHD="CHD_955_genes"
infix_DCM="DCM_909_genes"

data_subset_id="."
gene_subset_subset_list="."
gene_subset_subset_id="."

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

for infile in "${indir_gridss}"/*.gridss.*_genes.vpot_final_output_file.txt; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  family_and_samples="${array[0]}"
  IFS='_' read -r -a array <<< "$family_and_samples"
  family_or_sample="${array[0]}" # this is either the sample id or it is the family id if there is one
  len=${#array[@]}
  family=''
  sample=''
  if [[ $len -eq 1 ]]; then
    is_singleton=1
    is_family=0
    sample=$family_or_sample # =$family_and_samples
    string_for_grep="\t${sample}\t"
  else
    is_family=1
    is_singleton=0
    family=$family_or_sample
    string_for_grep="^${family}\t"
  fi

  # Do not process multi-sample files. They will have been turned into single-sample files. Cnvnator does not have multi-sample files.
  if [[ $sample != '' ]]; then

    gridss_infile=$(ls -1 "${indir_gridss}"/"${sample}".gridss.*_genes.vpot_final_output_file.txt | head -n 1)
    manta_infile=$(ls -1 "${indir_manta}"/"${sample}".manta.*_genes.vpot_final_output_file.txt | head -n 1)
    cnvnator_infile=$(ls -1 "${indir_cnvnator}"/"${sample}".cnvnator.*_genes.vpot_final_output_file.txt | head -n 1)

    # Use files that include large SVs > 1,000,000 bp size.
    # If a large variant is called by both cnvnator and gridss/manta, then it might be a true positive.

    echo -e "${sample}\t${outdir}\t${manta_infile}\t${gridss_infile}\t${cnvnator_infile}\t${data_subset_id}\t${gene_subset_subset_list}\t${gene_subset_subset_id}" >> "${out_manifest}"

  fi
done

echo 'out_manifest:' $out_manifest

./merge_results_from_structural_variants.sh $out_manifest

