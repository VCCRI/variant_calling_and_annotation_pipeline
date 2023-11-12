#!/bin/bash

# Run this file to create input commands to IGV
# Open IGV software ()
# In the top menu, navigate to Tools, then to Run Batch Scripts
# Open the file produced by this script
# The commands will run in IGV and automatically produce screenshots of the variants in the reports

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

bamdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
bamfile_suffix='.markDup.setTags.bqsr.bam'
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/igv_images
outfile_basename=all.commands_to_generate_igv_screenshots.annovar_clinvarDATE.txt
cohorts=(CHD DCM)
description='SNP_or_indel'

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

outfile="${outdir}"/"${outfile_basename}"
:>$outfile

for cohort in "${cohorts[@]}"; do

  echo 'Processing:' $cohort
  infile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples_extracts/"${cohort}"_"${batch}"__gatk_platypus_mnp.AllChrom.annovar_vep.possiblyPathogenic.vpot_final_output_file.vpot_gt_0.txt

  sample_col=-1
  gene_col=-1
  AN_col=-1
  FR_col=-1
  cols_string=$(head -n 1 $infile)
  IFS=$'\t' read -r -a array <<< "$cols_string"
  i=0
  for this_col in "${array[@]}"; do
    i=$(( i + 1 ))
    if [[ $this_col == "Sample_Id" ]]; then
      sample_col=$i
    fi
    if [[ $this_col == "Gene.refGene" ]]; then
      gene_col=$i
    fi
    if [[ $this_col == "INFO_Gene.refGene" ]]; then
      gene_col=$i
    fi
    if [[ $this_col == "AN" ]]; then
      AN_col=$i
    fi
    if [[ $this_col == "FR" ]]; then
      FR_col=$i
    fi
  done

  tmp_infile=tmp_infile_for_igv.txt
  echo 'grep -v' '^Rank' $infile '| awk -v sample_col='"$sample_col" '-v gene_col='"$gene_col" 'BEGIN {FS="\t";OFS="\t"} {print $sample_col, $3, $4, $5, $gene_col}' '| sort | uniq >' $tmp_infile '|| true'
  grep -v '^Rank' $infile | awk -v sample_col="$sample_col" -v gene_col="$gene_col" -v AN_col="$AN_col" -v FR_col="$FR_col" 'BEGIN {FS="\t";OFS="\t"} {
if ((($AN_col!=".")&&($AN_col<=1)) || (($FR_col!=".")&&($FR_col<0.1))){print $sample_col, $3, $4, $5, $gene_col}}' | sort | uniq > $tmp_infile || true
  current_sample=""

  while read sample chrom pos end gene; do

    if [[ "$sample" != "$current_sample" ]]; then
      current_sample=$sample
      outdir2="${outdir}"/"${sample}"
      mkdir -p "${outdir2}"
      inbam="${bamdir}"/"${sample}${bamfile_suffix}"

      echo '// batch file for IGV snapshot' >> $outfile
      echo 'new' >> $outfile
      echo 'genome' $genome_version '>>' $outfile
      echo 'setSleepInterval 30' >> $outfile
      echo 'maxPanelHeight 30000' >> $outfile
      echo 'snapshotDirectory' $outdir2 >> $outfile
      echo '//' >> $outfile
      echo 'new' >> $outfile
      echo 'load' $inbam >> $outfile
      echo '//' >> $outfile
    fi

    pos1=$(( pos - 1 ))
    pos2=$(( end + 1 ))
    coords="${chrom}":"${pos1}"-"${pos2}"
    outpng="${sample}"."${description}"."${chrom}"_"${pos}"."${gene}".png

    echo 'goto' $coords >> $outfile
    echo 'sort position' >> $outfile
    echo 'expand' >> $outfile
    echo 'snapshot' $outpng >> $outfile
    echo '//' >> $outfile

  done < $tmp_infile
  echo ''
done

echo ''
echo 'outfile:' $outfile
echo ''

