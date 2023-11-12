#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/platypus

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest="${outdir}"/"${cohort}".platypus.bam_input_list.txt

:>"${out_manifest}"
echo 'out_manifest:' $out_manifest

ls -1 "${indir}"/*.bam > $out_manifest

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

    jobid=plt"${chr_chrom}"

    outfile="${outdir}"/"${cohort}".platypus."${chr_chrom}".vcf

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

        echo 'qsub -N' ${jobid} '-v infile='"${out_manifest}"',outdir='"${outdir}"',outfile='"${outfile}"',chrom='"${chr_chrom}"',sw_and_refs='"${sw_and_refs}" 'callPlatypus_1_jointcalled_vcf.sh'
        qsub -N ${jobid} -v infile="${out_manifest}",outdir="${outdir}",outfile="${outfile}",chrom="${chr_chrom}",sw_and_refs="${sw_and_refs}" callPlatypus_1_jointcalled_vcf.sh && \
        touch "${queue_file}" && \
        echo "${jobid} queued" $out_manifest && \
        sleep 1

      fi
    #fi

done

