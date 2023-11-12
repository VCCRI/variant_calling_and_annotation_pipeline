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
while read sample infile outdir; do

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

  run_directory="${outdir}"

  jobid=spl"${sample}"

  outfile="${outdir}"/"${sample}".spliceogen_output

  queue_file="${outfile}.queued"
  lock_file="${outfile}.lock"
  done_file="${outfile}.done"
  term_file="${outfile}.term"
  log_file="${outfile}.log"

  # Spliceogen must run from a directory containing the spliceogena and maxentgen software that is being run.
  # So make a copy of the software in the outdir, if it is not already present.
  echo 'The Spliceogen software needs to be in the output directory in order to run Spliceogen.'
  echo 'Have you already run the following in order to copy Spliceogen software to the output directory?'
  echo ''
  echo 'cp -r' "${spliceogen_sw}"'/Spliceogen/*' "${outdir}"
  echo ''
  echo 'cp' "${spliceogen_sw}"'/resources/spliceogen_run_when_launched_in_pbs.sh' "${outdir}"
  echo ''

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

      qsub -N ${jobid} -v sample="${sample}",run_directory="${run_directory}",input="${infile}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}",gtf_for_spliceogen="${gtf_for_spliceogen}",fasta_for_spliceogen="${fasta_for_spliceogen}" spliceogen_1_launch_one_sample_to_pbs.sh # && \
      touch "${queue_file}" && \
      echo "${jobid} queued" $infile && \
      sleep 1
      num_queued=$((num_queued+1))

      echo ''
      echo ''
      echo ''
    fi
  #fi

done < $manifest
