#!/bin/bash

# All the samples in this manifest belong to the same family and will be joint-called together to produce one output file.
# The family id and outdir are taken from the first line of the manifest input file.

manifest=$1

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

infile_list=''
sample_list=''
outfile=''
jobid=''
out_family=''
out_outdir=''
family_and_samples=''

while read sample family infile outdir; do

  if [[ $infile_list == '' ]]; then
    infile_list=$infile
  else
    infile_list2="${infile_list}:${infile}"
    infile_list=$infile_list2
  fi

  if [[ $sample_list == '' ]]; then
    sample_list=$sample
  else
    sample_list2="${sample_list}:${sample}"
    sample_list=$sample_list2
  fi

  if [[ $out_family == '' ]]; then
    out_family=$family
  fi

  if [[ $family_and_samples == '' ]]; then
    family_and_samples="${family}_${sample}"
  else
    family_and_samples2="${family_and_samples}_${sample}"
    family_and_samples=$family_and_samples2
  fi

  if [[ $out_outdir == '' ]]; then
    out_outdir=$outdir
  fi

  outfile="${outdir}"/"${family_and_samples}".gridss.vcf.gz
  jobid=grd"${family_and_samples}"

done < $manifest

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

    echo 'qsub -N' ${jobid} '-v family='"${out_family}"',infile_list='"${infile_list}"',sample_list='"${sample_list}"',outdir='"${out_outdir}"',outfile='"${outfile}"',sw_and_refs='"${sw_and_refs}" 'gridss_1_jointcallSVs_gridss.sh'
    qsub -N ${jobid} -v family="${out_family}",infile_list="${infile_list}",sample_list="${sample_list}",outdir="${out_outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}" gridss_1_jointcallSVs_gridss.sh && \
    touch "${queue_file}" && \
    ##echo "${jobid} queued" $infile && \
    sleep 1
    num_queued=$((num_queued+1))
    echo ''

  fi
#fi
