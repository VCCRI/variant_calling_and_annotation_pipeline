#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/platypus
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/platypus

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

mito="M"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
fi

for chrom in {1..22} "X" "Y" "$mito" ; do

    chr_chrom="chr"$chrom
    if [[ $genome_version == "hg19" ]]; then
      chr_chrom=$chrom
    fi

    jobid=spl"${chrom}"

    infile="${indir}"/"${cohort}".platypus."${chr_chrom}".vcf
    outfile="${outdir}"/"${cohort}".platypus.vt_decompose."${chr_chrom}".vcf

    queue_file="${outfile}.queued"
    lock_file="${outfile}.lock"
    done_file="${outfile}.done"
    term_file="${outfile}.term"
    log_file="${outfile}.log"

    #if [ -e "${prev_done_file}" ]; then
      if [ -e "${queue_file}" ]; then
        echo "${jobid} already queued"
      elif [ -e "${lock_file}" ]; then
        echo "${jobid} already running" $lock_file
      elif [ -e "${done_file}" ]; then
        echo "${jobid} already done" $done_file
      elif [ -e "${term_file}" ]; then
        echo "${jobid} was terminated"
      else

        echo 'qsub -N' ${jobid} '-v infile='"${infile}"',outdir='"${outdir}"',outfile='"${outfile}"',chrom='"${chr_chrom}"',sw_and_refs='"${sw_and_refs}" 'callPlatypus_2_decompose.sh'
        qsub -N ${jobid} -v infile="${infile}",outdir="${outdir}",outfile="${outfile}",chrom="${chr_chrom}",sw_and_refs="${sw_and_refs}" callPlatypus_2_decompose.sh && \
        touch "${queue_file}" && \
        echo "${jobid} queued" $chr_chrom && \
        sleep 1

      fi
    #fi

done

