#!/bin/bash

manifest=$1

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

queue=normal
max_queued=$limit_num_of_pbs_submissions
num_queued=$((max_queued+1))

while read sample outdir in_manta in_gridss in_cnvnator data_subset_id gene_subset_subset_list gene_subset_subset_id; do

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

  outfile_prefix="${outdir}"/"${sample}".manta_gridss_cnvnator.merge
  if [[ "$data_subset_id" != "." ]]; then
    outfile_prefix="${outdir}"/"${sample}".manta_gridss_cnvnator.merge."${data_subset_id}"
  fi

  jobid=mrg"${sample}"

  queue_file="${outfile}.queued"
  lock_file="${outfile}.lock"
  done_file="${outfile}.done"
  term_file="${outfile}.term"
  log_file="${outfile}.log"

  #if [ -e "${prev_done_file}" ]; then
  #  if [ -e "${queue_file}" ]; then
  #    echo "${jobid} already queued"
  #  elif [ -e "${lock_file}" ]; then
  #    echo "${jobid} already running"
  #  elif [ -e "${done_file}" ]; then
  #    echo "${jobid} already done"
  #  elif [ -e "${term_file}" ]; then
  #    echo "${jobid} was terminated"
  #  else

    module load R/3.6.1
    module unload intel-fc intel-cc
    module load intel-compiler/2019.3.199
    

    echo 'Rscript' "${sw}"'/victorchang_scripts/merge_results_from_structural_variants.R' $outfile_prefix $in_manta $in_gridss $in_cnvnator $gene_subset_subset_list $gene_subset_subset_id
    Rscript "${sw}"/victorchang_scripts/merge_results_from_structural_variants.R $outfile_prefix $in_manta $in_gridss $in_cnvnator $gene_subset_subset_list $gene_subset_subset_id

    echo ''
    echo 'Sizes of output files:'
    wc -l "${outfile_prefix}"*
    echo ''; echo ''; echo ''

  #  fi
  #fi

done < $manifest
