#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

maf=0.01
infile_prefix=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate/"${batch}".gatk_vqsr_and_platypus_mnp.chr
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.vcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_vpot_ranking/
# MY_BATCH.gatk_vqsr_and_platypus_mnp.chr1.annovar_clinvarDATE.hg38_multianno.vep_reformatted.extras.all_coding_regions.maf_0p01.vcf

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

param_file="${currdir}"/default_VPOT_MAF_0p01_PPF_for_CHD_names.txt

mito22Y="22MY"
if [[ $genome_version == "hg19" ]]; then
  mito22Y="22MTY"
fi


# 22/M/Y have been merged to ensure that all samples have at least one variant in the file, otherwise vpot drops samples and variants
for chrom in {1..21} "X" "$mito22Y" ; do

    infile="${infile_prefix}${chrom}${infile_suffix}"
    infile_basename=$(basename $infile)

    maf_as_text="${maf/\./p}"
    outfile_basename="${infile_basename%.vcf}".maf_"${maf_as_text}"_and_scoring.vcf
    outfile="${outdir}"/"${outfile_basename}"

    jobid=maf"${chrom}"

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

        echo 'qsub -N' ${jobid} '-v param_file='"${param_file}"',infile='"${infile}"',outdir='"${outdir}"',outfile='"${outfile}"',sw_and_refs='"${sw_and_refs}" 'vpotInheritance_2_filter_maf_chrom_for_vpot.sh'
        qsub -N ${jobid} -v param_file="${param_file}",infile="${infile}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}" vpotInheritance_2_filter_maf_chrom_for_vpot.sh # && \
        #touch "${queue_file}" && \
        #echo "${jobid} queued" $infile && \
        #sleep 1

        #touch "${queue_file}"
        #echo './vpotInheritance_2_filter_maf_chrom_for_vpot.sh' $param_file $infile $outdir $outfile $sw_and_refs
        #./vpotInheritance_2_filter_maf_chrom_for_vpot.sh $param_file $infile $outdir $outfile $sw_and_refs

    #  fi
    #fi
done

