#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l mem=2G
#PBS -l ncpus=1
#PBS -l jobfs=2G
#PBS -N SVextract
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile_tsv=$2
#outdir=$3
#outfile=$4
#outprefix=$5
#sw_and_refs=$6
#cohort=$7
#currdir=$8
#large_homozygous_runs=$9
#homozygous_runs="${10}"

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

touch "${lock_file}"
rm -f "${queue_file}"

module load python3
module load bedtools
module load R/3.6.1
module unload intel-fc intel-cc
module load intel-compiler/2019.3.199


# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

tmpout_exon_hits_tsv="${tmpdir}"/tmpout_exon_hits_tsv.tsv

##### Identify exon columns and extract all records where they are not empty.

col_gene_CDSexons=0
col_gencode_gene_CDSexons=0

echo 'col_gene_CDSexons_result=head -n 1' $infile_tsv '| sed -e s/\t/\n/g | grep -n gene_CDSexons'
col_gene_CDSexons_result=`head -n 1 $infile_tsv | sed -e 's/\t/\n/g' | grep -n gene_CDSexons`
echo ''
echo 'col_gencode_gene_CDSexons_result=head -n 1' $infile_tsv '| sed -e s/\t/\n/g | grep -n gencode_gene_CDSexons || true'
col_gencode_gene_CDSexons_result=`head -n 1 $infile_tsv | sed -e 's/\t/\n/g' | grep -n gencode_gene_CDSexons || true`
echo ''

IFS=':' read -r -a array <<< "$col_gene_CDSexons_result"
col_gene_CDSexons="${array[0]}"

if [[ $col_gencode_gene_CDSexons_result == "" ]]; then
  col_gencode_gene_CDSexons=$col_gene_CDSexons
else
  IFS=':' read -r -a array <<< "$col_gencode_gene_CDSexons_result"
  col_gencode_gene_CDSexons="${array[0]}"
fi

echo 'col_gene_CDSexons =' $col_gene_CDSexons
echo 'col_gencode_gene_CDSexons =' $col_gencode_gene_CDSexons

echo 'awk -v col_gene_CDSexons='$col_gene_CDSexons '-v col_gencode_gene_CDSexons='$col_gencode_gene_CDSexons 'BEGIN {FS="\t";OFS="\t"} {if ( (($col_gene_CDSexons!=".")&&($col_gene_CDSexons!="")) || (($col_gencode_gene_CDSexons!=".")&&($col_gencode_gene_CDSexons!="")) ) {print$0}}' $infile_tsv '>' $tmpout_exon_hits_tsv
awk -v col_gene_CDSexons=$col_gene_CDSexons -v col_gencode_gene_CDSexons=$col_gencode_gene_CDSexons 'BEGIN {FS="\t";OFS="\t"} {if ( (($col_gene_CDSexons!=".")&&($col_gene_CDSexons!="")) || (($col_gencode_gene_CDSexons!=".")&&($col_gencode_gene_CDSexons!="")) ) {print$0}}' $infile_tsv > $tmpout_exon_hits_tsv
echo ''


##### Add sample id and gene flags and incidentalome flags

temp_outfile_0="${tmpdir}"/temp_outfile_0.txt
temp_outfile_1="${tmpdir}"/temp_outfile_1.txt
temp_outfile_2="${tmpdir}"/temp_outfile_2.txt
temp_outfile_3="${tmpdir}"/temp_outfile_3.txt
temp_outfile_4="${tmpdir}"/temp_outfile_4.txt

echo ''
echo 'Add sample id.'
echo 'awk -v sample='"$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmpout_exon_hits_tsv}" '>' "${temp_outfile_0}"

awk -v sample="$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmpout_exon_hits_tsv}" > "${temp_outfile_0}"
echo ''

report_gene_column_names='gene_CDSexons,gene_entireGene,FANTOM5_heart_TSS_gene,gencode_gene_CDSexons,gencode_gene_entireGene'

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

temp_outfile_incidentalome_1="${tmpdir}"/temp_outfile_incidentalome_1.txt
temp_outfile_incidentalome_2="${tmpdir}"/temp_outfile_incidentalome_2.txt

temp_outfile_large_hrun_flag="${tmpdir}"/temp_outfile_large_hrun_flag.txt
temp_outfile_hrun_flag="${tmpdir}"/temp_outfile_hrun_flag.txt

echo ''
echo 'Add VCCRI incidentalome.'
echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $tmp_outfile_with_gene_flags $temp_outfile_incidentalome_1 $name1 $genes1 '.' $report_gene_column_names
echo ''

Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $tmp_outfile_with_gene_flags $temp_outfile_incidentalome_1 $name1 $genes1 . $report_gene_column_names

echo ''
echo 'Add VCGS incidentalome.'
echo 'Rscript' $sw'/victorchang_scripts/add_gene_flag_to_report.R' $temp_outfile_incidentalome_1 $temp_outfile_incidentalome_2 $name2 $genes2 '.' $report_gene_column_names
echo ''

Rscript $sw/victorchang_scripts/add_gene_flag_to_report.R $temp_outfile_incidentalome_1 $temp_outfile_incidentalome_2 $name2 $genes2 . $report_gene_column_names

echo 'awk -f' $sw'/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file='"$large_homozygous_runs" '-v annot='"proband_BigHomozygRun" $temp_outfile_incidentalome_2 '>' $temp_outfile_large_hrun_flag
awk -f $sw/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file="$large_homozygous_runs" -v annot="proband_BigHomozygRun" $temp_outfile_incidentalome_2 > $temp_outfile_large_hrun_flag
echo ''

echo 'awk -f' $sw'/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file='"$homozygous_runs" '-v annot='"proband_HRun" $temp_outfile_large_hrun_flag '>' $temp_outfile_hrun_flag
awk -f $sw/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file="$homozygous_runs" -v annot="proband_HRun" $temp_outfile_large_hrun_flag > $temp_outfile_hrun_flag
echo ''

echo ''
echo 'cp' $temp_outfile_hrun_flag $outfile
cp $temp_outfile_hrun_flag $outfile

echo ''
echo 'outfile:' $outfile
echo ''

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"



