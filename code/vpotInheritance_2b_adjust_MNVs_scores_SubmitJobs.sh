#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

maf=0.01

infile_prefix=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_vpot_ranking/"${batch}".gatk_vqsr_and_platypus_mnp.chr
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.maf_0p01_and_scoring.vcf_final_output_file.txt

outfile_prefix=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_vpot_ranking/"${batch}".gatk_vqsr_and_platypus_mnp.chr
outfile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.maf_0p01_and_scoring.vcf_final_output_file.adjust_scores.txt

outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_vpot_ranking

mito22Y="22MY"
if [[ $genome_version == "hg19" ]]; then
  mito22Y="22MTY"
fi

# 22/M/Y have been merged to ensure that all samples have at least one variant in the file, otherwise vpot drops samples and variants
for chrom in {1..21} "X" "$mito22Y" ; do

    infile="${infile_prefix}${chrom}${infile_suffix}"
    outfile="${outfile_prefix}${chrom}${outfile_suffix}"

    jobid=adj"${chrom}"

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

        echo 'qsub -N' ${jobid} '-v infile='"${infile}"',outdir='"${outdir}"',outfile='"${outfile}"',sw_and_refs='"${sw_and_refs}" 'vpotInheritance_2b_adjust_MNVs_scores.pbs'
        qsub -N ${jobid} -v infile="${infile}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}" vpotInheritance_2b_adjust_MNVs_scores.pbs # && \
        #touch "${queue_file}" && \
        #echo "${jobid} queued" $infile && \
        #sleep 1
        echo ''

    #  fi
    #fi
done

