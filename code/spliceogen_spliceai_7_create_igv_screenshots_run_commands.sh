#!/bin/bash

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/spliceogen_spliceai
bamdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
bamfile_suffix=.markDup.setTags.bqsr.bam

outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/igv_images_for_splicing
outfile="${outdir}"/all.commands_to_generate_igv_screenshots.splicing.txt
:>$outfile

cohorts=(CHD DCM)
#cohorts=(DCM)

for cohort in "${cohorts[@]}"; do

 echo 'Processing:' $cohort

 for infile in "${indir}"/MY_BATCH__merge_annovar*."${cohort}"_*.tsv; do

  echo 'infile' $infile

  sample_col=-1
  gene_col=-1
  AN_col=-1
  FR_col=-1
  spliceogen_donGainP=-1
  spliceogen_accGainP=-1
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
    if [[ $this_col == "spliceogen_donGainP" ]]; then
      spliceogen_donGainP_col=$i
    fi
    if [[ $this_col == "spliceogen_accGainP" ]]; then
      spliceogen_accGainP_col=$i
    fi
  done
  echo 'spliceogen_donGainP_col='$spliceogen_donGainP_col
  echo 'spliceogen_accGainP_col='$spliceogen_accGainP_col

  # sort output by spliceogen_donGainP and spliceogen_accGainP in case there are too many to do all of them

  tmp_infile=tmp_infile_for_igv.txt
  echo 'grep -v' '^Rank' $infile '| awk -v sample_col='"$sample_col" '-v gene_col='"$gene_col" '-v spliceogen_donGainP_col='"$spliceogen_donGainP_col" '-v spliceogen_accGainP_col='"$spliceogen_accGainP_col" 'BEGIN {FS="\t";OFS="\t"} {print $spliceogen_donGainP_col, $spliceogen_accGainP_col, $sample_col, $1, $2, $3, $gene_col}' '| sort | uniq | cut -d$''\t' '-f3- >' $tmp_infile '|| true'
  grep -v '^Rank' $infile | awk -v sample_col="$sample_col" -v gene_col="$gene_col" -v spliceogen_donGainP_col="$spliceogen_donGainP_col" -v spliceogen_accGainP_col="$spliceogen_accGainP_col" -v AN_col="$AN_col" -v FR_col="$FR_col" 'BEGIN {FS="\t";OFS="\t"} {if ((($AN_col!="")&&($AN_col<=1)) || (($FR_col!="")&&($FR_col<0.1))) {print $spliceogen_donGainP_col, $spliceogen_accGainP_col, $sample_col, $1, $2, $3, $gene_col}}' | sort | uniq | cut -d$'\t' -f3- > $tmp_infile || true
  current_sample=""

  while read sample chrom pos end gene; do

    if [[ "$sample" != "$current_sample" ]]; then
      current_sample=$sample
      outdir2="${outdir}"/"${sample}"
      mkdir -p "${outdir2}"
      inbam="${bamdir}"/"${sample}${bamfile_suffix}"

      echo '// batch file for IGV snapshot' >> $outfile
      echo 'new' >> $outfile
      echo 'genome hg38' >> $outfile
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
    outpng="${sample}"."${chrom}"_"${pos}"."${gene}".png

    echo 'goto' $coords >> $outfile
    echo 'sort position' >> $outfile
    echo 'expand' >> $outfile
    echo 'snapshot' $outpng >> $outfile
    echo '//' >> $outfile

  done < $tmp_infile
  echo ''
 done
done

echo ''
echo 'outfile:' $outfile
echo ''


