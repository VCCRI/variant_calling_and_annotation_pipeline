#!/bin/bash

manifest=$1

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

queue=normal
max_queued=$limit_num_of_pbs_submissions
num_queued=$((max_queued+1))

while read sample infile_tsv outdir; do

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

  outfile="${outdir}"/"${sample}".cnvnator.filter.tsv

  jobid=fil"${sample}"

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

        #qsub -N ${jobid} -v sample="${sample}",infile_tsv="${infile_tsv}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}" cnvnator_2_filter.sh # && \
        #touch "${queue_file}" && \
        #echo "${jobid} queued" $infile && \
        #sleep 1
        #num_queued=$((num_queued+1))

	./cnvnator_2_filter.sh $sample $infile_tsv $outdir $outfile $sw_and_refs

      fi
  #fi

done < $manifest

