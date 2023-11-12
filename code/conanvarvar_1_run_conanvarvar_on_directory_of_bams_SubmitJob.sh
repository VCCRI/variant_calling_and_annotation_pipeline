#!/bin/bash

manifest=$1

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

queue=normal
max_queued=250
num_queued=$((max_queued+1))

in_sw_dir=$conanvarvar_sv_dir

# This job runs ConanVarvar to process an entire directory of bam files in one job.
# The ConanVarvar software needs to be copied and run from is own directory because the software directory is also a working directory.

while read in_bam_dir outdir; do

  # mkdir -p "${outdir}"

  # Only submit 100 jobs at a time
  #if [ ${num_queued} -gt ${max_queued} ]; then
  #  set +e
  #  num_queued=$(qstat -u `whoami` | grep "${queue}" | awk '($10 == "Q"){total+=1} END{print total+0}')
  #  set -e
  #fi
  #if [ ${num_queued} -gt ${max_queued} ]; then
  #  echo -e "Terminating: queue full"
  #  break 10
  #fi

  outdir_sw_and_working_dir="${outdir}"/conanvarvar_sw_and_working_dir
  outdir_results="${outdir}"/conanvarvar_results

  mkdir -p "${outdir_sw_and_working_dir}"
  mkdir -p "${outdir_results}"

  cp -r "${in_sw_dir}"/* "${outdir_sw_and_working_dir}"

  jobid=conanvarvar

  qsub -N ${jobid} -v in_bam_dir="${in_bam_dir}",outdir_sw_and_working_dir="${outdir_sw_and_working_dir}",outdir_results="${outdir_results}",sw_and_refs="${sw_and_refs}" conanvarvar_1_run_conanvarvar_on_directory_of_bams.sh

done < $manifest
