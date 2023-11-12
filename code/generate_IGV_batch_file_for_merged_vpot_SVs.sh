#!/bin/bash
set -euo pipefail

indir1=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate
indir2=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
indir3=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate
inbamdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/igv_images_for_SVs
outfile="${outdir}"/all_samples_SV.generate_IGV_batch_file_for_SVs.txt

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

infile_list=$(ls -1 "${indir1}"/*_"${batch}"__*_genes.vpot_final_output_file*.txt "${indir2}"/*_"${batch}"__*_genes.vpot_final_output_file*.txt "${indir3}"/*_"${batch}"__*_genes.vpot_final_output_file*.txt)

tmpdir=./tmp_for_generate_IGV_batch_file_for_merged_vpot_SVs
mkdir -p $tmpdir
tmpfile1="${tmpdir}"/generate_IGV_batch_file_for_vpot_SVs.tmpfile1.txt
tmpfile2="${tmpdir}"/generate_IGV_batch_file_for_vpot_SVs.tmpfile2.txt

# Make a list of variants and their sample (and family). Write to a temporary file. It is not ordered by sample.

# CHD_family_MY_BATCH__gridss.CHD_955_genes.vpot_final_output_file.txt
# CHD_singles_MY_BATCH__gridss.CHD_955_genes.vpot_final_output_file.vpot_gt_0.txt
# DCM_MY_BATCH__gridss.DCM_909_genes.vpot_final_output_file.vpot_gt_0.txt
:>$tmpfile1
for infile in $infile_list; do

  echo 'infile' $infile
  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  cohort_subset_tool="${array[0]}"
  #echo 'cohort_subset_tool' $cohort_subset_tool

  cols_string=$(head -n 1 $infile)
  IFS=$'\t' read -r -a array <<< "$cols_string"
  i=0
  col_chrom=-1
  col_pos=-1
  col_end=-1
  col_svtype=-1
  col_gene=-1
  col_sample=-1
  col_family_sample_ids=-1
  col_match1=-1
  for this_col in "${array[@]}"; do
    i=$(( i + 1 ))
    if [[ $this_col == "#CHROM" ]]; then
      col_chrom=$i
    fi
    if [[ $this_col == "POS" ]]; then
      col_pos=$i
    fi
    if [[ $this_col == "END" ]]; then
      col_end=$i
    fi
    if [[ $this_col == "SVTYPE" ]]; then
      col_svtype=$i
    fi
    if [[ $this_col == "cnv_type" ]]; then
      col_svtype=$i
    fi
    if [[ $this_col == "gene_entireGene" ]]; then
      col_gene=$i
    fi
    if [[ $this_col == "Sample_Id" ]]; then
      col_sample=$i
    fi
    if [[ $this_col == "family_sample_ids" ]]; then
      col_family_sample_ids=$i
    fi
    if [[ $this_col == "otherSamples_90percentMatch" ]]; then
      col_match1=$i
    fi
    if [[ $this_col == "includeLowQual_otherSamples_90percentMatch" ]]; then
      col_match1=$i
    fi
  done

  echo 'col_match1' $col_match1
  #echo 'col_family_sample_ids' $col_family_sample_ids
  #echo 'col_sample' $col_sample
  if [[ "$col_family_sample_ids" -eq -1 ]]; then
    awk -v col_family_sample_ids="$col_family_sample_ids" -v col_sample="$col_sample" -v col_chrom="$col_chrom" -v col_pos="$col_pos" -v col_end="$col_end" -v col_svtype="$col_svtype" -v col_gene="$col_gene" -v col_match1="$col_match1" 'BEGIN {FS="\t";OFS="\t"} {if (($1!="Ranking") && (($5-$4)<500000) && ($1>=0.75) && ($col_svtype!="BND") && ($col_match1==".")) {print ".", $col_sample, $col_chrom, $col_pos, $col_end, $col_svtype, $col_gene}}' $infile >> $tmpfile1
  else
    awk -v col_family_sample_ids="$col_family_sample_ids" -v col_sample="$col_sample" -v col_chrom="$col_chrom" -v col_pos="$col_pos" -v col_end="$col_end" -v col_svtype="$col_svtype" -v col_gene="$col_gene" -v col_match1="$col_match1" 'BEGIN {FS="\t";OFS="\t"} {if (($1!="Ranking") && (($5-$4)<500000) && ($1>=0.75) && ($col_svtype!="BND") && ($col_match1==".")) {print $col_family_sample_ids, $col_sample, $col_chrom, $col_pos, $col_end, $col_svtype, $col_gene}}' $infile >> $tmpfile1
  fi
  #echo ''
done

# Sort the temporary file of list of variants by family and sample, and remove variants because there are too many.

sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 $tmpfile1 > $tmpfile2
echo $tmpfile2
wc -l $tmpfile2

# Write out the IGV commands to view these variants.
# At the beginning of each new family or sample, write out IGV commands to open new bam files.

IFS=$'\t'
:>$outfile
prev_family_or_sample=""
while read family_sample_ids sample chrom pos end svtype gene; do

  #echo 'family_sample_ids' $family_sample_ids
  #echo 'sample' $sample
  #echo 'chrom' $chrom
  #echo 'pos' $pos
  #echo 'end' $end
  #echo 'svtype' $svtype
  #echo 'gene' $gene
  process_this_var=1
  if [[ $end == "END" ]]; then
    process_this_var=0
  else
    svlen=$(( $end - $pos ))
    if [[ $svlen -gt 1000000 ]]; then
      process_this_var=0
    fi
  fi

  if [[ $process_this_var == "1" ]]; then

    IFS=',' read -r -a array <<< "$gene"
    gene="${array[0]}"
    IFS=$'\t'
    if [[ $gene == "." ]]; then
      gene=""
    fi

    if [[ $family_sample_ids == "." ]]; then
      curr_family_or_sample=$sample
    else
      curr_family_or_sample=$family_sample_ids
   fi

    # Write commands for IGV to load new bam(s)
    if [[ $curr_family_or_sample != $prev_family_or_sample ]]; then

      if [[ $family_sample_ids == "." ]]; then
        # for just one bam
        outdir2="${outdir}"/$sample
        #mkdir -p $outdir2
        inbam="${inbamdir}"/"${sample}".markDup.setTags.bqsr.bam
        echo 'inbam' $inbam
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
      else
        # for multiple bams, because this is a family/trio
        outdir2="${outdir}"/$family_sample_ids
        #mkdir -p $outdir2
        echo '// batch file for IGV snapshot' >> $outfile
        echo 'new' >> $outfile
        echo 'genome hg38' >> $outfile
        echo 'setSleepInterval 30' >> $outfile
        echo 'maxPanelHeight 30000' >> $outfile
        echo 'snapshotDirectory' $outdir2 >> $outfile
        echo '//' >> $outfile
        echo 'new' >> $outfile
        IFS='_' read -r -a array <<< "$family_sample_ids"
        IFS=$'\t'
        for element in "${array[@]}"; do
          char3=${element:0:3}
          if [[ $char3 != "FAM" ]]; then
            inbam="${inbamdir}"/"${element}".markDup.setTags.bqsr.bam
            echo 'load' $inbam >> $outfile
          fi
        done
        echo '//' >> $outfile
      fi
    fi

    # Write commands for IGV to go to coords and print out a png
    pos_minus_100=$(( $pos - 100 ))
    pos_plus_100=$(( $pos + 100 ))
    coords="${chrom}:${pos_minus_100}-${pos_plus_100}"
    outpng="${curr_family_or_sample}_${chrom}_${pos_minus_100}_${pos_plus_100}_${svtype}_${gene}".png
    echo 'goto' $coords >> $outfile
    echo 'sort position' >> $outfile
    echo 'expand' >> $outfile
    echo 'snapshot' $outpng >> $outfile
    echo '//' >> $outfile
    pos_minus_100=$(( $end - 100 ))
    pos_plus_100=$(( $end + 100 ))
    coords="${chrom}:${pos_minus_100}-${pos_plus_100}"
    outpng="${curr_family_or_sample}_${chrom}_${pos_minus_100}_${pos_plus_100}_${svtype}_${gene}".png
    echo 'goto' $coords >> $outfile
    echo 'sort position' >> $outfile
    echo 'expand' >> $outfile
    echo 'snapshot' $outpng >> $outfile
    echo '//' >> $outfile
  fi

  if [[ $family_sample_ids == "." ]]; then
    prev_family_or_sample=$sample
  else
    prev_family_or_sample=$family_sample_ids
  fi
  #echo ''
done < $tmpfile2

echo ''
echo 'output:' $outfile
echo ''



