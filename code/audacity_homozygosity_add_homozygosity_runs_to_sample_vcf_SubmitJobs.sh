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
while read sample infile outdir flagA large_homozygous_runs_file flagB homozygous_runs_file sample_job_differentiating_id; do

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
  outfile_basename="${infile_basename%.gz}"
  outfile_basename="${outfile_basename%.vcf}"
  outfile_basename="${outfile_basename%.hg38_multianno}"
  outfile_basename="${outfile_basename%.hg19_multianno}"
  outfile_basename="${outfile_basename}".homrun.vcf

  outfile="${outdir}"/"${outfile_basename}"

  jobid=vhm"${sample}"_"${sample_job_differentiating_id}"

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

      echo 'qsub -N' ${jobid} '-v sample='"${sample}"',infile_vcf='"${infile}"',infile2_flag='"${flagA}"',infile2='"${large_homozygous_runs_file}"',infile3_flag='"${flagB}"',infile3='"${homozygous_runs_file}"',outdir='"${outdir}"',outfile_basename='"${outfile_basename}"',sw_and_refs='"${sw_and_refs}"',sample_job_differentiating_id='"${sample_job_differentiating_id}" 'audacity_homozygosity_add_homozygosity_runs_to_sample_vcf.sh'
      qsub -N ${jobid} -v sample="${sample}",infile_vcf="${infile}",infile2_flag="${flagA}",infile2="${large_homozygous_runs_file}",infile3_flag="${flagB}",infile3="${homozygous_runs_file}",outdir="${outdir}",outfile_basename="${outfile_basename}",sw_and_refs="${sw_and_refs}",sample_job_differentiating_id="${sample_job_differentiating_id}" audacity_homozygosity_add_homozygosity_runs_to_sample_vcf.sh # && \
      touch "${queue_file}" && \
      echo "${jobid} queued" $infile_tsv # && \
      sleep 1
      num_queued=$((num_queued+1))

      #echo './audacity_homozygosity_add_homozygosity_runs_to_sample_vcf.sh' $sample $infile $flagA $large_homozygous_runs_file $flagB $homozygous_runs_file $outdir $outfile_basename $sw_and_refs $sample_job_differentiating_id
      #./audacity_homozygosity_add_homozygosity_runs_to_sample_vcf.sh $sample $infile $flagA $large_homozygous_runs_file $flagB $homozygous_runs_file $outdir $outfile_basename $sw_and_refs $sample_job_differentiating_id

      echo ''
    fi
  #fi

done < $manifest
