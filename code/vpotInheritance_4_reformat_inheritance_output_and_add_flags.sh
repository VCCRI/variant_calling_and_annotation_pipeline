#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -l jobfs=10MB
#PBS -N vpotInheritance
#PBS -lstorage=gdata/abcd

set -euo pipefail

#infile=$1
#outdir=$2
#outfile=$3
#sw_and_refs=$4
#cohort=$5

echo ''
echo 'infile:' $infile
echo 'outdir:' $outdir
echo 'outfile:' $outfile
echo 'sw_and_refs:' $sw_and_refs
echo ''

out_prefix=$outfile

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

term_handler()
{
    rm -f "${lock_file}"
    touch "${term_file}"
    exit 1
}
trap 'term_handler' TERM

touch "${lock_file}"
rm -f "${queue_file}"

module unload intel-fc intel-cc
module load intel-mkl/2019.3.199
module load python3
module load R/3.6.1
module load intel-compiler/2019.3.199


# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

currdir=$(dirname $sw_and_refs)

echo 'gawk -f' $sw'/victorchang_scripts/convert_file_with_semicolon_info_to_tab_delimited_include_all_info_headers.awk' $infile '>' $outfile
gawk -f $sw/victorchang_scripts/convert_file_with_semicolon_info_to_tab_delimited_include_all_info_headers.awk $infile > $outfile
echo ''

  temp_prefix="${tmpdir}"/temp_prefix
  temp_outfile_1="${temp_prefix}".temp_outfile_1.txt
  temp_outfile_2="${temp_prefix}".temp_outfile_2.txt
  temp_outfile_3="${temp_prefix}".temp_outfile_3.txt
  temp_outfile_4="${temp_prefix}".temp_outfile_4.txt

  report_gene_column_names='Gene.refGene,Gene.'"$wgEncodeGencodeBasic"',spliceai_gene,spliceaiIndels_gene,FANTOM5_CAGEr_humanHeart_transcriptStartSites'

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

  tmp_outfile_with_gene_flags=$outfile

  if [[ $name1 != "" ]]; then
    echo ''
    echo 'Add first flag.' $cohort $name1
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_1 $name1 $genes1 $regions1 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_1 $name1 $genes1 $regions1 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_1
  fi

  if [[ $name2 != "" ]]; then
    echo ''
    echo 'Add second flag.' $cohort $name2
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_2 $name2 $genes2 $regions2 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_2 $name2 $genes2 $regions2 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_2
  fi

  if [[ $name3 != "" ]]; then
    echo ''
    echo 'Add third flag.' $cohort $name3
    echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_3 $name3 $genes3 $regions3 $report_gene_column_names
    echo ''
    Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_3 $name3 $genes3 $regions3 $report_gene_column_names
    tmp_outfile_with_gene_flags=$temp_outfile_3
  fi

  if [[ $name4 != "" ]]; then
    echo ''
    echo 'Add fourth flag.' $cohort $name4
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

  temp_outfile_large_hrun_flag="${temp_prefix}".temp_outfile_large_hrun_flag.txt
  temp_outfile_hrun_flag="${temp_prefix}".temp_outfile_hrun_flag.txt

  outfile_flags="${outfile%.tsv}"
  outfile_flags="${outfile_flags%.txt}"
  outfile_flags="${outfile_flags%.tab_delimited}".flags.tsv

  echo ''
  echo 'Add VCCRI incidentalome.'
  echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_incidentalome_1 $name1 $genes1 '.' $report_gene_column_names
  echo ''

  Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_incidentalome_1 $name1 $genes1 . $report_gene_column_names

  echo ''
  echo 'Add VCGS incidentalome.'
  echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $temp_outfile_incidentalome_1 $outfile_flags $name2 $genes2 '.' $report_gene_column_names
  echo ''

  Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $temp_outfile_incidentalome_1 $outfile_flags $name2 $genes2 . $report_gene_column_names

echo ''
echo 'outfile:' $outfile
echo ''
echo 'outfile_flags:' $outfile_flags
echo ''
echo 'Finished!'

touch "${done_file}"
rm -f "${lock_file}"

