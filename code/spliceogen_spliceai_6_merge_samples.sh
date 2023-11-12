#!/bin/bash

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

cohorts=(CHD DCM)

for cohort in "${cohorts[@]}"; do

  echo 'Processing:' $cohort

  temp_hdr="${indir}"/temp_hdr_"${cohort}".txt
  extract_desc=''
  if [[ $cohort == "CHD" ]]; then
    extract_desc='CHD_955_genes'
  fi
  if [[ $cohort == "DCM" ]]; then
    extract_desc='DCM_909_genes'
  fi

  ########## /my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai/SAMPLE_ONE.annovar_clinvarDATE_spliceai.spliceogen.gnomad_filtered_splicing_variants.CHD_955_genes.tsv

  outfile="${indir}"/"${batch}"__merge_annovar_clinvar_spliceai.vep.spliceogen.gnomad_filtered_splicing_variants."${extract_desc}".tsv
  all_files="${indir}"/*.annovar_clinvar*_spliceai.vep.spliceogen.gnomad_filtered_splicing_variants."${extract_desc}".tsv

  highest_count=0
  for infile in $all_files; do
    count=`head -n 1 $infile | sed -e 's/\t/\n/g' | wc -l | cut -d' ' -f1`
    if [[ "$count" -gt "$highest_count" ]]; then
      highest_count=$count
      one_file=$infile
    fi
  done
  head -n 1 $one_file > $temp_hdr

  col1=0
  col2=0
  echo 'col1_result=head -n 1' $one_file '| sed -e s/\t/\n/g | grep -P -n "sai_site|spliceai_site" | head -n 1'
  col1_result=`head -n 1 $one_file | sed -e 's/\t/\n/g' | grep -P -n "sai_site|spliceai_site" | head -n 1`|| true
  echo ''
  echo 'col2_result=head -n '1 $one_file '| sed -e s/\t/\n/g | grep -P -n "saiI_site|spliceaiIndels_site" | head -n 1'
  col2_result=`head -n 1 $one_file | sed -e 's/\t/\n/g' | grep -P -n "saiI_site|spliceaiIndels_site" | head -n 1`|| true
  echo ''
  IFS=':' read -r -a array <<< "$col1_result"
  col1="${array[0]}"
  IFS=':' read -r -a array <<< "$col2_result"
  col2="${array[0]}"
  echo 'col1 =' $col1
  echo 'col2 =' $col2
  echo ''

  echo 'cat' "${all_files}" '| grep -v ^chrom | grep -v ^CHROM | LC_ALL=C sort -t$\t -k'"${col1}"','"${col1}"'r -k'"${col2}"','"${col2}"'r | cat' $temp_hdr '- >' $outfile
  cat $all_files | grep -v '^chrom' | grep -v '^CHROM' | LC_ALL=C sort -t$'\t' -k"${col1}","${col1}"r -k"${col2}","${col2}"r | cat $temp_hdr - > $outfile || true
  echo ''

  echo 'outfile:' $outfile
  echo ''

  ########## /my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai/SAMPLE_ONE.annovar_clinvarDATE_spliceai.spliceogen.gnomad_filtered_splicing_variants.CHD_955_genes.needs_spliceai_tensorflow_scores.tensorflow_results.tsv

  outfile="${indir}"/"${batch}"__merge_annovar_clinvar.vep.spliceogen.gnomad_filtered_splicing_variants."${extract_desc}".merge_spliceai_tensorflow_scores.tsv
  all_files="${indir}"/*.annovar_clinvar*.vep.spliceogen.gnomad_filtered_splicing_variants."${extract_desc}".merge_spliceai_tensorflow_scores.tsv

  highest_count=0
  for infile in $all_files; do
    count=`head -n 1 $infile | sed -e 's/\t/\n/g' | wc -l | cut -d' ' -f1`
    if [[ "$count" -gt "$highest_count" ]]; then
      highest_count=$count
      one_file=$infile
    fi
  done
  head -n 1 $one_file > $temp_hdr

  col1=0
  col2=0
  echo 'col1_result=head -n 1' $one_file '| sed -e s/\t/\n/g | grep -P -n "SpliceAI" | head -n 1'
  col1_result=`head -n 1 $one_file | sed -e 's/\t/\n/g' | grep -P -n "SpliceAI" | head -n 1`|| true
  echo ''
  IFS=':' read -r -a array <<< "$col1_result"
  col1="${array[0]}"
  echo 'col1 =' $col1
  echo ''

  echo 'cat' "${all_files}" '| grep -v ^chrom | grep -v ^CHROM | LC_ALL=C sort -t$\t -k'"${col1}"','"${col1}"'r | cat' $temp_hdr '- >' $outfile
  cat $all_files | grep -v '^chrom' | grep -v '^CHROM' | LC_ALL=C sort -t$'\t' -k"${col1}","${col1}"r | cat $temp_hdr - > $outfile || true
  echo ''

  echo '##################################################'
  echo 'outfile:' $outfile
  echo '##################################################'
  echo ''

done

echo ''
echo 'Finished!'
echo ''

