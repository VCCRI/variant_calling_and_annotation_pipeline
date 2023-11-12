#!/bin/bash

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vpot_inheritance
# FAMILY_FOUR.maf_0p01.AllChrom.vpot_inheritance.AD.AD_father_SAMPLE_ONE_AD_variant_filtered_output_file.txt
# FAMILY_FOUR.maf_0p01.AllChrom.vpot_inheritance.AD.AD_mother_SAMPLE_ONE_AD_variant_filtered_output_file.txt
# FAMILY_FOUR.maf_0p01.AllChrom.vpot_inheritance.AR.AR_SAMPLE_ONE_AR_variant_filtered_output_file.txt
# FAMILY_FOUR.maf_0p01.AllChrom.vpot_inheritance.CH.CH_SAMPLE_ONE_CH_variant_filtered_output_file.txt
# FAMILY_FOUR.maf_0p01.AllChrom.vpot_inheritance.DN.DN_SAMPLE_ONE_DN_variant_filtered_output_file.txt
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vpot_inheritance

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

pedigree_file="${batch}".pedigree_file.ped
in_manifest="${batch}".one_line_per_family.txt

for infile in "${indir}"/*.vpot_inheritance.*_variant_filtered_output_file.txt; do # covers CH DN AR AD CaseControl

    infile_basename=$(basename $infile)
    outfile="${infile%.txt}".tab_delimited.tsv

    IFS='.' read -r -a array <<< "$infile_basename"
    family="${array[0]}"
    IFS='.' read -r -a array <<< "$family"
    family="${array[0]}"

    echo 'family:' $family
    family_for_grep="\t${family}\t|\t${family}$"
    samples_cohort=`grep -P "${family_for_grep}" "${in_manifest}"`
    IFS=$'\t' read -r -a array <<< "$samples_cohort"
    cohort="${array[0]}"

    jobid=inf"${infile_basename}"

    queue_file="${outfile}.queued"
    lock_file="${outfile}.lock"
    done_file="${outfile}.done"
    term_file="${outfile}.term"
    log_file="${outfile}.log"

    #if [ -e "${prev_done_file}" ]; then

      if [ -e "${queue_file}" ]; then
        echo "${jobid} already queued" $queue_file
      elif [ -e "${lock_file}" ]; then
        echo "${jobid} already running" $lock_file
      elif [ -e "${done_file}" ]; then
        echo "${jobid} already done" $done_file
      elif [ -e "${term_file}" ]; then
        echo "${jobid} was terminated" $term_file
      else

          echo 'qsub -N' ${jobid} '-v infile='"${infile}"',outdir='"${outdir}"',outfile='"${outfile}"',sw_and_refs='"${sw_and_refs}"',cohort='"${cohort}" 'vpotInheritance_4_reformat_inheritance_output_and_add_flags.sh'
          qsub -N ${jobid} -v infile="${infile}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}",cohort="${cohort}" vpotInheritance_4_reformat_inheritance_output_and_add_flags.sh # && \
          touch "${queue_file}" && \
          echo "${jobid} queued" $infile && \
          sleep 1

          echo ''
      fi
    #fi

done

