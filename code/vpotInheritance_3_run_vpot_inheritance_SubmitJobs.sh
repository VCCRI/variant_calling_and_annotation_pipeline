#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

maf=0.01

infile_prefix=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_vpot_ranking/"${batch}".gatk_vqsr_and_platypus_mnp.
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.maf_0p01_and_scoring.vcf_final_output_file.adjust_scores.txt

outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vpot_inheritance

maf_as_text="${maf/\./p}"

pedigree_file="${currdir}"/"${batch}".pedigree_file.ped

IFS=$'\n' list_families=($(awk 'BEGIN {FS="\t";OFS="\t"} {if ( (($3!=0) || ($4!=0)) && ($1!="Family") ) {print $0}}' $pedigree_file | cut -d$'\t' -f1 | sort | uniq))

chrom="AllChrom"

# AD = autosomal dominant
# AR = autosomal recessive
# CH = compound hets
# DN = denovo

for family in "${list_families[@]}"; do
    echo 'family' $family

    for inheritance_model in "AD" "AR" "CH" "DN" "CaseControl"; do

      infile="${infile_prefix}""${chrom}${infile_suffix}"
      outfile_basename="${family}".maf_"${maf_as_text}"."${chrom}".vpot_inheritance."${inheritance_model}"
      outfile="${outdir}"/"${outfile_basename}"

      jobid=inh"${chrom}${inheritance_model}"_"${family}"

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

          echo 'qsub -N' ${jobid} '-v family='"${family}"',inheritance_model='"${inheritance_model}"',infile='"${infile}"',pedigree_file='"${pedigree_file}"',maf='"${maf}"',outdir='"${outdir}"',outfile='"${outfile}"',sw_and_refs='"${sw_and_refs}" 'vpotInheritance_3_run_vpot_inheritance.pbs'
          qsub -N ${jobid} -v family="${family}",inheritance_model="${inheritance_model}",infile="${infile}",pedigree_file="${pedigree_file}",maf="${maf}",outdir="${outdir}",outfile="${outfile}",sw_and_refs="${sw_and_refs}" vpotInheritance_3_run_vpot_inheritance.pbs # && \
          touch "${queue_file}"
          echo "${jobid} queued" $infile
          sleep 1
          echo ''

          #touch "${queue_file}"
          #echo './vpotInheritance_3_run_vpot_inheritance_nopbs.sh' $family $inheritance_model $infile $pedigree_file $maf $outdir $outfile $sw_and_refs
          #./vpotInheritance_3_run_vpot_inheritance_nopbs.sh $family $inheritance_model $infile $pedigree_file $maf $outdir $outfile $sw_and_refs

          echo ''
        fi
      #fi
    done
done


