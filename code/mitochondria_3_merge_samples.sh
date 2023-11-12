#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

manifest="${batch}".one_line_per_family.txt
indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/mitochondria
infile_basename_suffix=.chromM.helixmtdb.mitimpact.mitotip.mitomap.extramito.possibly_pathogenic.tsv
outfile_basename_suffix=_"${batch}"__mitochondria.possibly_pathogenic.tsv
tmphdr='temp_hdr.txt'
tmphdr2='temp_hdr2.txt'

cohorts=(CHD DCM)

for cohort in "${cohorts[@]}"; do
  echo 'Processing:' $cohort
  outfile="${indir}"/"${cohort}${outfile_basename_suffix}"
  string_of_files=''

  for infile in "${indir}"/*"${infile_basename_suffix}"; do

    infile_basename=$(basename $infile)
    IFS='.' read -r -a array <<< "$infile_basename"
    this_sample="${array[0]}"
    this_family="${array[0]}"

    sample_for_grep="\t${this_sample}\t|\t${this_sample}$"
    samples_cohort=`grep -P "${sample_for_grep}" "${manifest}"`
    IFS=$'\t' read -r -a array <<< "$samples_cohort"
    this_cohort="${array[0]}"

    if [[ "$this_cohort" == "$cohort" ]]; then
      echo 'Add' $infile
      new_string_of_files="${string_of_files} ${infile}"
      string_of_files=$new_string_of_files
      one_file=$infile
    fi
  done

  head -n 1 $one_file > $tmphdr
  last_col=$(sed -e 's/\t/\n/g' "$tmphdr" | tail -n 1)
  sed_string='s/'"${last_col}"'$/sample_field/'
  sed $sed_string $tmphdr > $tmphdr2
  cat $string_of_files | grep -v '^CHROM' | sort -k1nr -k2nr | cat $tmphdr2 - > $outfile
  echo ''
  echo 'outfile:' $outfile
  echo ''
done
