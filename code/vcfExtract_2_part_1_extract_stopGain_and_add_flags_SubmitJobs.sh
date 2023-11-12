#!/bin/bash

manifest=$1

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

queue=normal
max_queued=$limit_num_of_pbs_submissions
num_queued=$((max_queued+1))

# one loop per input file in the the input manifest file
while read sample infile_prefix infile_suffix outdir cohort; do # only genes required by cohort will be output

  mkdir -p "${outdir}"

  if [[ $limit_num_of_pbs_submissions -gt 0 ]]; then
    # Only submit 100 jobs at a time
    if [ ${num_queued} -gt ${max_queued} ]; then
      set +e
      num_queued=$(qstat -u `whoami` | grep "${queue}" | awk '($10 == "Q") {total+=1} END {print total+0}')
      set -e
    fi
    if [ ${num_queued} -gt ${max_queued} ]; then
      echo -e "Terminating: queue full"
      break 10
    fi
  fi

  infile_prefix_basename=$(basename $infile_prefix)
  infile_prefix_basename="${infile_prefix_basename%.chr}"
  infile_suffix_prefix="${infile_suffix%.gz}"
  infile_suffix_prefix="${infile_suffix_prefix%.vcf}"

  outfile="${outdir}"/"${infile_prefix_basename}".AllChrom"${infile_suffix_prefix}".stoploss_stopgain.tsv

  jobid=stp"${sample}"

  queue_file="${outfile}.queued"
  lock_file="${outfile}.lock"
  done_file="${outfile}.done"
  term_file="${outfile}.term"
  log_file="${outfile}.log"

  #if [ -e "${prev_done_file}" ]; then
    if [ -e "${queue_file}" ]; then
      echo "${jobid} already queued"
    elif [ -e "${lock_file}" ]; then
      echo "${jobid} already running"
    elif [ -e "${done_file}" ]; then
      echo "${jobid} already done"
    elif [ -e "${term_file}" ]; then
      echo "${jobid} was terminated"
    else

      qsub -N ${jobid} -v sample="${sample}",infile_prefix="${infile_prefix}",infile_suffix="${infile_suffix}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}",cohort="${cohort}" vcfExtract_2_part_1_extract_stopGain_and_add_flags.sh # && \
      touch "${queue_file}" && \
      echo "${jobid} queued" $infile # && \
      sleep 1
      num_queued=$((num_queued+1))

      #./vcfExtract_2_part_1_extract_stopGain_and_add_flags.sh $sample $infile_prefix $infile_suffix $outdir $outfile $sw_and_refs $cohort

    fi
  #fi

done < $manifest


