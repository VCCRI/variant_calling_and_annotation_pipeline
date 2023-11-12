#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -l jobfs=16G
#PBS -N extracts
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile_prefix=$2
#infile_suffix=$3
#outdir=$4
#outfile=$5
#sw_and_refs=$6
#cohort=$7

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

touch "${lock_file}"
rm -f "${queue_file}"

module load R/3.6.1
module unload intel-fc intel-cc
module load intel-compiler/2019.3.199


# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

random_number=$RANDOM
script_name=$(basename $0)
infile_prefix_basename=$(basename $infile_prefix)
infile_suffix_prefix="${infile_suffix%.gz}"
infile_suffix_prefix="${infile_suffix_prefix%.vcf}"
outfile_basename=$(basename $outfile)
outfile_basename_prefix="${outfile_basename%.tsv}"

tmpdir="${PBS_JOBFS}"/tmp/vcfExtract_2_part_1."${outfile_basename}"
mkdir -p "${tmpdir}"

currdir=$(dirname $sw_and_refs)

rm -rf "${outfile}"

echo ''
echo 'Extract stoploss, stopgain and frameshift variants from whole genome for sample' $infile_prefix'*'$infile_suffix
echo ''

tmp_outfile_prev_chroms="${tmpdir}"/tmp_outfile_prev_chroms.tsv

for chrom in {1..22} "M" "X" "Y" ; do

  echo ''
  echo 'Processing chrom' $chrom
  echo ''

  infile_chrom_vcf="${infile_prefix}${chrom}${infile_suffix}"
  infile_chrom_vcf_basename=$(basename $infile_chrom_vcf)
  infile_chrom_vcf_basename_prefix="${infile_chrom_vcf_basename%.gz}"
  infile_chrom_vcf_basename_prefix="${infile_chrom_vcf_basename_prefix%.vcf}"

  tmpout_chrom_stoploss_stopgain_vcf="${tmpdir}"/"${infile_chrom_vcf_basename_prefix}".stoploss_stopgain_frameshift.vcf
  tmpout_chrom_stoploss_stopgain_tsv="${tmpdir}"/"${infile_chrom_vcf_basename_prefix}".stoploss_stopgain_frameshift.tsv

  echo ''
  echo 'Extract stoploss/stopgain/frameshift variants from' $infile_chrom_vcf 'to produce' $tmpout_chrom_stoploss_stopgain_vcf 'and after that' $tmpout_chrom_stoploss_stopgain_tsv
  echo ''

  echo ''
  echo 'grep -P -i ^#|transcript_ablation|splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|stoploss|stopgain|\tframeshift' $infile_chrom_vcf '>' $tmpout_chrom_stoploss_stopgain_vcf
  echo ''

  set +e
  set +o pipefail
  grep -P -i '^#|transcript_ablation|splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|stop_lost|start_lost|transcript_amplification|stoploss|stopgain|\tframeshift' $infile_chrom_vcf > $tmpout_chrom_stoploss_stopgain_vcf
  set -e
  set -e pipefail

  echo ''
  echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i' "${tmpout_chrom_stoploss_stopgain_vcf}" '-o' "${tmpout_chrom_stoploss_stopgain_tsv}" '-end_id YES'
  echo ''

  python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i "${tmpout_chrom_stoploss_stopgain_vcf}" -o "${tmpout_chrom_stoploss_stopgain_tsv}" -end_id YES

  ##### Add sample id and gene flags and incidentalome flags

  infile_tsv_basename=$(basename $tmpout_chrom_stoploss_stopgain_tsv)
  temp_prefix="${tmpdir}"/temp_prefix
  temp_outfile_0="${temp_prefix}".temp_outfile_0.txt
  temp_outfile_1="${temp_prefix}".temp_outfile_1.txt
  temp_outfile_2="${temp_prefix}".temp_outfile_2.txt
  temp_outfile_3="${temp_prefix}".temp_outfile_3.txt
  temp_outfile_4="${temp_prefix}".temp_outfile_4.txt

  echo ''
  echo 'Chrom:' $chrom 'Add sample id.'
  echo 'awk -v sample='"$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmpout_chrom_stoploss_stopgain_tsv}" '>' "${temp_outfile_0}"
  echo ''

  awk -v sample="$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmpout_chrom_stoploss_stopgain_tsv}" > "${temp_outfile_0}"

  report_gene_column_names='INFO_Gene.refGene,INFO_Gene.'"$wgEncodeGencodeBasic"',INFO_spliceai_gene,spliceaiIndels_gene,INFO_FANTOM5_CAGEr_humanHeart_transcriptStartSites'

  name1=""
  name2=""
  name3=""
  name4=""

  if [[ "$cohort" == "CHD" ]]; then

    regions1="${currdir}"/list_genes_CHD_tier1_106_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes1="${currdir}"/list_genes_CHD_tier1_106_genes_withoutBars.txt
    name1=$(basename $genes1)
    name1="${name1%.txt}"
    name1="${name1%_withoutBars}"
    name1="${name1/list_genes_/}"

    regions2="${currdir}"/list_genes_CHD_tier2_227_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes2="${currdir}"/list_genes_CHD_tier2_227_genes_withoutBars.txt
    name2=$(basename $genes2)
    name2="${name2%.txt}"
    name2="${name2%_withoutBars}"
    name2="${name2/list_genes_/}"

    regions3="${currdir}"/list_genes_CHD_eCHD_403_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes3="${currdir}"/list_genes_CHD_eCHD_403_genes_withoutBars.txt
    name3=$(basename $genes3)
    name3="${name3%.txt}"
    name3="${name3%_withoutBars}"
    name3="${name3/list_genes_/}"

    regions4="${currdir}"/list_genes_CHD_nonCHDtier1tier2_207_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes4="${currdir}"/list_genes_CHD_nonCHDtier1tier2_207_genes_withoutBars.txt
    name4=$(basename $genes4)
    name4="${name4%.txt}"
    name4="${name4%_withoutBars}"
    name4="${name4/list_genes_/}"

  fi

  if [[ "$cohort" == "DCM" ]]; then

    regions1="${currdir}"/list_genes_DCM_tier1_list_of_32_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes1="${currdir}"/list_genes_DCM_tier1_list_of_32_genes_withoutBars.txt
    name1=$(basename $genes1)
    name1="${name1%.txt}"
    name1="${name1%_withoutBars}"
    name1="${name1/list_genes_/}"

    regions2="${currdir}"/list_genes_DCM_tier2_list_of_14_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes2="${currdir}"/list_genes_DCM_tier2_list_of_14_genes_withoutBars.txt
    name2=$(basename $genes2)
    name2="${name2%.txt}"
    name2="${name2%_withoutBars}"
    name2="${name2/list_genes_/}"

    regions3="${currdir}"/list_genes_DCM_other_Flagship_nonTier1nonTier2_list_of_527_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes3="${currdir}"/list_genes_DCM_other_Flagship_nonTier1nonTier2_list_of_527_genes_withoutBars.txt
    name3=$(basename $genes3)
    name3="${name3%.txt}"
    name3="${name3%_withoutBars}"
    name3="${name3/list_genes_/}"

    regions4="${currdir}"/list_genes_DCM_Minoche2018_and_extra_516_genes_hg38_RefSeq_regions_tab_delimited.txt
    genes4="${currdir}"/list_genes_DCM_Minoche2018_and_extra_516_genes_withoutBars.txt
    name4=$(basename $genes4)
    name4="${name4%.txt}"
    name4="${name4%_withoutBars}"
    name4="${name4/list_genes_/}"

  fi

  tmp_outfile_with_gene_flags=$temp_outfile_0

  if [[ $name1 != "" ]]; then
    echo ''
    echo 'Chrom:' $chrom 'Add first flag.' $cohort $name1
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_1 $name1 $genes1 $regions1 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_1 $name1 $genes1 $regions1 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_1
  fi

  if [[ $name2 != "" ]]; then
    echo ''
    echo 'Chrom:' $chrom 'Add second flag.' $cohort $name2
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_2 $name2 $genes2 $regions2 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_2 $name2 $genes2 $regions2 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_2
  fi

  if [[ $name3 != "" ]]; then
    echo ''
    echo 'Chrom:' $chrom 'Add third flag.' $cohort $name3
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_3 $name3 $genes3 $regions3 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_3 $name3 $genes3 $regions3 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_3
  fi

  if [[ $name4 != "" ]]; then
    echo ''
    echo 'Chrom:' $chrom 'Add fourth flag.' $cohort $name4
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_4 $name4 $genes4 $regions4 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_4 $name4 $genes4 $regions4 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_4
  fi

  name1='VCCRI_Incidentalome'
  genes1="${currdir}"/list_genes_VCCRI_Incidentalome_withoutBars.txt
  name2='VCGS_Incidentalome'
  genes2="${currdir}"/list_genes_VCGS_Incidentalome_genes_withoutBars.txt

  temp_outfile_incidentalome_1="${temp_prefix}".temp_outfile_incidentalome_1.txt
  temp_outfile_incidentalome_2="${temp_prefix}".temp_outfile_incidentalome_2.txt

  echo ''
  echo 'Chrom:' $chrom 'Add VCCRI incidentalome.'
  echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_incidentalome_1 $name1 $genes1 '.' $report_gene_column_names
  echo ''

  Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_incidentalome_1 $name1 $genes1 . $report_gene_column_names

  echo ''
  echo 'Chrom:' $chrom 'Add VCGS incidentalome.'
  echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $temp_outfile_incidentalome_1 $temp_outfile_incidentalome_2 $name2 $genes2 '.' $report_gene_column_names
  echo ''

  Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $temp_outfile_incidentalome_1 $temp_outfile_incidentalome_2 $name2 $genes2 . $report_gene_column_names

  ##### Concatenate this chromosome output to the output of previous chromosomes.

  if [[ -f "${outfile}" ]]; then

    echo ''
    echo 'Chrom:' $chrom 'Make a copy of the stopGain/stopLoss/frameshift of previous chromosomes.'
    echo 'cp' $outfile $tmp_outfile_prev_chroms
    echo ''
    cp $outfile $tmp_outfile_prev_chroms

    echo ''
    echo 'Chrom:' $chrom 'Concatenate the stopGain/stopLoss/frameshift variants in previous chromosomes to the current chromosome just created.'
    echo 'grep -v "^CHROM"' $temp_outfile_incidentalome_2 '| cat' $tmp_outfile_prev_chroms '- >' $outfile
    echo ''
    grep -v "^CHROM" $temp_outfile_incidentalome_2 | cat $tmp_outfile_prev_chroms - > $outfile || true

  else

    echo ''
    echo 'Chrom:' $chrom 'This is the first chromosome so this results of stopGain/stopLoss/frameshift variants in only matching genes is the currently running output.'
    echo 'cp' $temp_outfile_incidentalome_2 $outfile
    echo ''
    cp $temp_outfile_incidentalome_2 $outfile

  fi

done

echo ''
echo 'outfile:' $outfile
echo ''

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

