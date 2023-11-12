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
while read sample infile infile2 outdir; do

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

  infile_basename=$(basename $infile)
  infile_basename_prefix="${infile_basename%.txt}"
  infile_basename_prefix="${infile_basename_prefix%.tsv}"
  infile_basename_prefix="${infile_basename_prefix%.tensorflow_results}"
  infile_basename_prefix="${infile_basename_prefix%.needs_spliceai_tensorflow_scores}"

  outfile="${outdir}"/"${infile_basename_prefix}".merge_spliceai_tensorflow_scores.tsv

  jobid=smg"${chrom}"

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

      #echo 'qsub -N' ${jobid} '-v sample='"${sample}"',infile='"${infile}"',infile2='"${infile2}"',outdir='"${outdir}"',outfile='"${outfile}"',sw_and_refs='"${sw_and_refs}" 'spliceogen_spliceai_4b_add_spliceai_tensorflow_output_to_tsv.sh'
      #qsub -N ${jobid} -v sample="${sample}",infile="${infile}",infile2="${infile2}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}" spliceogen_spliceai_4b_add_spliceai_tensorflow_output_to_tsv.sh && \
      #touch "${queue_file}" && \
      #echo "${jobid} queued" $infile && \
      #sleep 1
      #num_queued=$((num_queued+1))

      echo './spliceogen_spliceai_4b_add_spliceai_tensorflow_output_to_tsv.sh' $sample $infile $infile2 $outdir $outfile $sw_and_refs
      ./spliceogen_spliceai_4b_add_spliceai_tensorflow_output_to_tsv.sh $sample $infile $infile2 $outdir $outfile $sw_and_refs

      echo ''
    fi
  #fi

done < $manifest