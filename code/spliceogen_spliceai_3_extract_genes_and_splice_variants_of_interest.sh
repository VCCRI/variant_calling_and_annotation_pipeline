#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=64G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -N spliceaiextract
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile=$2
#outdir=$3
#outfile=$4
#sw_and_refs=$5
#cohort=$6

module load R/3.6.1

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

outfile_basename=$(basename $outfile)

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

genes=""
regions=""
currdir=$(dirname $sw_and_refs)

if [[ "$cohort" == "CHD" ]]; then
  regions="${currdir}"/list_genes_CHD_955_genes_hg38_RefSeq_regions_colon_dash_format.txt
  genes="${currdir}"/list_genes_CHD_955_genes_withoutBars.txt
  outfile="${outfile%.tsv}"
  outfile="${outfile%.txt}"
  outfile="${outfile%.extracts}"
  outfile="${outfile%.collapsed}"
  outfile="${outfile}".CHD_955_genes.tsv
fi

if [[ "$cohort" == "DCM" ]]; then
  regions="${currdir}"/list_genes_DCM_909_genes_hg38_RefSeq_regions_colon_dash_format.txt
  genes="${currdir}"/list_genes_DCM_909_genes_withoutBars.txt
  outfile="${outfile%.tsv}"
  outfile="${outfile%.txt}"
  outfile="${outfile%.extracts}"
  outfile="${outfile%.collapsed}"
  outfile="${outfile}".DCM_909_genes.tsv
fi

outfile2="${outfile%.tsv}".needs_spliceai_tensorflow_scores.tsv

##### List all the genes for a first rough pass.

echo ''
echo 'Extract variants from' $infile 'having input genes' $genes
echo ''

# Selecting by gene name can bring back too many hits, so select by gene name in context of tsv file.
# 
#	Gene.refGene=ABCB5				INFO_Gene.refGene
#	Gene.refGene==LINC02564\x3bLOC102723376
#	Gene.wgEncodeGencodeBasicV33=ABCB5		INFO_Gene.wgEncodeGencodeBasicV33
#	Gene.wgEncodeGencodeBasicV33=NONE\x3bAC215217.1
#	Gene.wgEncodeGencodeBasicV33=TUBB8\x3bZMYND11
#	GENE=ABCB5
#	CLINVAR20200629_GENEINFO=BRCA1:672
#	CLINVAR20200629_GENEINFO=CASQ2:845|VANGL1:81839
#	CLINVAR20200629_GENEINFO=TNNT2:7139
#	sai_gene=ABCB5
#	FANTOM5_CAGEr_humanHeart_transcriptStartSites=Name\x3dX237627570_237627644_+_240.4977261_LRRFIP1
#	FANTOM5_CAGEr_humanHeart_transcriptStartSites=Name\x3dX201365214_201365314_-_0_TNNT2,AC119427.2,X201365186_201365281_+_202.7365015_
#			AC119427.2 is an error in annotation reference file, needs to be changed to 201365186_201365281_+_202.7365015_AC119427.2

genes_string_for_filtering_tsv='^#|^CHROM|^chrom|^Chrom|^CHR\t|^Chr\t|^chr\t'
while IFS= read -r inline; do

  this_genes_string_for_filtering_tsv="\\t${inline}|;${inline}|\\x3b${inline}|\|${inline}|,${inline}|_${inline}|=${inline}"
  this_genes_string_for_filtering_tsv_2="${this_genes_string_for_filtering_tsv}|${inline}\\t|${inline};|${inline}\\x3b|${inline}:|${inline},|${inline}\|"

  new_genes_string_for_filtering_tsv="${genes_string_for_filtering_tsv}|${this_genes_string_for_filtering_tsv_2}"
  genes_string_for_filtering_tsv=$new_genes_string_for_filtering_tsv

done < "$genes"

tmpout_rough_infile_tsv="${tmpdir}"/tmpout_rough_infile_tsv.tsv
tmpout_genes_exact_tsv="${tmpdir}"/tmpout_genes_exact_tsv.tsv
tmpout_genes_exact_filtered_splicing_tsv="${tmpdir}"/tmpout_genes_exact_filtered_splicing_tsv.tsv
tmpout_genes_multinucleotide_splicing_tsv="${tmpdir}"/tmpout_genes_multinucleotide_splicing_tsv.tsv

##### Do a rough extract of fields that have a splicing entry, so that input is not so large
echo ''
echo 'Do a rough extract of fields that have a splicing entry, so that input is not so large'
echo ''

col_withinSite=0
col_donGainP=0
col_accGainP=0
col_donLossP=0
col_accLossP=0
spliceai_acceptor_gain=0
spliceai_acceptor_loss=0
spliceai_donor_gain=0
spliceai_donor_loss=0
spliceaiIndels_acceptor_gain=0
spliceaiIndels_acceptor_loss=0
spliceaiIndels_donor_gain=0
spliceaiIndels_donor_loss=0
mmsplice_delta_logit_psi=0

echo 'col_withinSite_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n spliceogen_withinSite'
col_withinSite_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n spliceogen_withinSite`
echo ''
echo 'col_donGainP_result=head -n '1 $infile '| sed -e s/\t/\n/g | grep -n spliceogen_donGainP'
col_donGainP_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n spliceogen_donGainP`
echo ''
echo 'col_accGainP_result=head -n '1 $infile '| sed -e s/\t/\n/g | grep -n spliceogen_accGainP'
col_accGainP_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n spliceogen_accGainP`
echo ''
echo 'col_donLossP_result=head -n '1 $infile '| sed -e s/\t/\n/g | grep -n spliceogen_donLossP'
col_donLossP_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n spliceogen_donLossP`
echo ''
echo 'col_accLossP_result=head -n '1 $infile '| sed -e s/\t/\n/g | grep -n spliceogen_accLossP'
col_accLossP_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n spliceogen_accLossP`
echo ''
echo 'col_spliceai_acceptor_gain_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n sai_AG_DS'
col_spliceai_acceptor_gain_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n sai_AG_DS`
echo ''
echo 'col_spliceai_acceptor_loss_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n sai_AL_DS'
col_spliceai_acceptor_loss_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n sai_AL_DS`
echo ''
echo 'col_spliceai_donor_gain_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n sai_DG_DS'
col_spliceai_donor_gain_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n sai_DG_DS`
echo ''
echo 'col_spliceai_donor_loss_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n sai_DL_DS'
col_spliceai_donor_loss_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n sai_DL_DS`
echo ''
echo 'col_spliceaiIndels_acceptor_gain_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n saiI_AG_DS'
col_spliceaiIndels_acceptor_gain_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n saiI_AG_DS`
echo ''
echo 'col_spliceaiIndels_acceptor_loss_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n saiI_AL_DS'
col_spliceaiIndels_acceptor_loss_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n saiI_AL_DS`
echo ''
echo 'col_spliceaiIndels_donor_gain_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n saiI_DG_DS'
col_spliceaiIndels_donor_gain_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n saiI_DG_DS`
echo ''
echo 'col_spliceaiIndels_donor_loss_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n saiI_DL_DS'
col_spliceaiIndels_donor_loss_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n saiI_DL_DS`
echo ''
echo 'col_mmsplice_delta_logit_psi_result=head -n 1' $infile '| sed -e s/\t/\n/g | grep -n mmsplice_delta_logit_psi'
col_mmsplice_delta_logit_psi_result=`head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n mmsplice_delta_logit_psi`
echo ''

IFS=':' read -r -a array <<< "$col_withinSite_result"
col_withinSite="${array[0]}"
IFS=':' read -r -a array <<< "$col_donGainP_result"
col_donGainP="${array[0]}"
IFS=':' read -r -a array <<< "$col_accGainP_result"
col_accGainP="${array[0]}"
IFS=':' read -r -a array <<< "$col_donLossP_result"
col_donLossP="${array[0]}"
IFS=':' read -r -a array <<< "$col_accLossP_result"
col_accLossP="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceai_acceptor_gain_result"
col_spliceai_acceptor_gain="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceai_acceptor_loss_result"
col_spliceai_acceptor_loss="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceai_donor_gain_result"
col_spliceai_donor_gain="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceai_donor_loss_result"
col_spliceai_donor_loss="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceaiIndels_acceptor_gain_result"
col_spliceaiIndels_acceptor_gain="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceaiIndels_acceptor_loss_result"
col_spliceaiIndels_acceptor_loss="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceaiIndels_donor_gain_result"
col_spliceaiIndels_donor_gain="${array[0]}"
IFS=':' read -r -a array <<< "$col_spliceaiIndels_donor_loss_result"
col_spliceaiIndels_donor_loss="${array[0]}"
IFS=':' read -r -a array <<< "$col_mmsplice_delta_logit_psi_result"
col_mmsplice_delta_logit_psi="${array[0]}"

echo 'col_withinSite =' $col_withinSite
echo 'col_donGainP =' $col_donGainP
echo 'col_accGainP =' $col_accGainP
echo 'col_donLossP =' $col_donLossP
echo 'col_accLossP =' $col_accLossP
echo 'col_spliceai_acceptor_gain =' $col_spliceai_acceptor_gain
echo 'col_spliceai_acceptor_loss =' $col_spliceai_acceptor_loss
echo 'col_spliceai_donor_gain =' $col_spliceai_donor_gain
echo 'col_spliceai_donor_loss =' $col_spliceai_donor_loss
echo 'col_spliceaiIndels_acceptor_gain =' $col_spliceaiIndels_acceptor_gain
echo 'col_spliceaiIndels_acceptor_loss =' $col_spliceaiIndels_acceptor_loss
echo 'col_spliceaiIndels_donor_gain =' $col_spliceaiIndels_donor_gain
echo 'col_spliceaiIndels_donor_loss =' $col_spliceaiIndels_donor_loss
echo 'col_mmsplice_delta_logit_psi =' $col_mmsplice_delta_logit_psi
echo ''

echo 'awk -v col_withinSite='$col_withinSite '-v donGainP='$col_donGainP '-v accGainP='$col_accGainP '-v donLossP' $col_donLossP '-v accLossP='$col_accLossP '-v col_spliceai_acceptor_gain='$col_spliceai_acceptor_gain '-v col_spliceai_acceptor_loss='$col_spliceai_acceptor_loss '-v col_spliceai_donor_gain='$col_spliceai_donor_gain '-v col_spliceai_donor_loss='$col_spliceai_donor_loss '-v col_spliceaiIndels_acceptor_gain='$col_spliceaiIndels_acceptor_gain '-v col_spliceaiIndels_acceptor_loss='$col_spliceaiIndels_acceptor_loss '-v col_spliceaiIndels_donor_gain='$col_spliceaiIndels_donor_gain '-v col_spliceaiIndels_donor_loss='$col_spliceaiIndels_donor_loss '-v col_mmsplice_delta_logit_psi='$col_mmsplice_delta_logit_psi 'BEGIN {FS="\t";OFS="\t"} {if ( (($col_withinSite!=".")&&(($col_donGainP!=".")||($col_accGainP!=".")||($col_donLossP!=".")||($col_accLossP!="."))) || ($col_spliceai_acceptor_gain!=".")||($col_spliceai_acceptor_loss!=".")||($col_spliceai_donor_gain!=".")||($col_spliceai_donor_loss!=".")||($col_spliceaiIndels_acceptor_gain!=".")||($col_spliceaiIndels_acceptor_loss!=".")||($col_spliceaiIndels_donor_gain!=".")||($col_spliceaiIndels_donor_loss!=".")||($col_mmsplice_delta_logit_psi!=".") ) {print $0}}' $infile '>' $tmpout_rough_infile_tsv
awk -v col_withinSite=$col_withinSite -v donGainP=$col_donGainP -v accGainP=$col_accGainP -v donLossP=$col_donLossP -v accLossP=$col_accLossP -v col_spliceai_acceptor_gain=$col_spliceai_acceptor_gain -v col_spliceai_acceptor_loss=$col_spliceai_acceptor_loss -v col_spliceai_donor_gain=$col_spliceai_donor_gain -v col_spliceai_donor_loss=$col_spliceai_donor_loss -v col_spliceaiIndels_acceptor_gain=$col_spliceaiIndels_acceptor_gain -v col_spliceaiIndels_acceptor_loss=$col_spliceaiIndels_acceptor_loss -v col_spliceaiIndels_donor_gain=$col_spliceaiIndels_donor_gain -v col_spliceaiIndels_donor_loss=$col_spliceaiIndels_donor_loss -v col_mmsplice_delta_logit_psi=$col_mmsplice_delta_logit_psi 'BEGIN {FS="\t";OFS="\t"} {if ( (($col_withinSite!=".")&&(($col_donGainP!=".")||($col_accGainP!=".")||($col_donLossP!=".")||($col_accLossP!="."))) || ($col_spliceai_acceptor_gain!=".")||($col_spliceai_acceptor_loss!=".")||($col_spliceai_donor_gain!=".")||($col_spliceai_donor_loss!=".")||($col_spliceaiIndels_acceptor_gain!=".")||($col_spliceaiIndels_acceptor_loss!=".")||($col_spliceaiIndels_donor_gain!=".")||($col_spliceaiIndels_donor_loss!=".")||($col_mmsplice_delta_logit_psi!=".") ) {print $0}}' $infile > $tmpout_rough_infile_tsv
echo ''

##### Extract variants in the genes of interest, using the infile as input
echo ''
echo 'Extract variants in the genes of interest, using the infile as input'
echo ''

echo ''
echo 'Rscript' "${sw}"'/victorchang_scripts/extract_by_optional_gene_name.R' "${tmpout_rough_infile_tsv}" "${tmpout_genes_exact_tsv}" "${genes}"
echo ''

Rscript "${sw}"/victorchang_scripts/extract_by_optional_gene_name.R "${tmpout_rough_infile_tsv}" "${tmpout_genes_exact_tsv}" "${genes}"

##### Extract splicing variants from spliceogen and/or spliceai, filter by gnomad, using the extracted genes of interest as input
echo ''
echo 'Extract splicing variants from spliceogen and/or spliceai, filter by gnomad, using the extracted genes of interest as input'
echo ''

filter_gnomad=0.05

echo ''
echo 'Rscript' "${sw}"'/victorchang_scripts/extract_spliceogen_and_spliceai_variants.R' "${tmpout_genes_exact_tsv}" "${tmpout_genes_exact_filtered_splicing_tsv}" "${filter_gnomad}"
echo ''

Rscript "${sw}"/victorchang_scripts/extract_spliceogen_and_spliceai_variants.R "${tmpout_genes_exact_tsv}" "${tmpout_genes_exact_filtered_splicing_tsv}" "${filter_gnomad}"

##### Add sample id and gene flags and incidentalome flags
echo ''
echo 'Add sample id and gene flags and incidentalome flags'
echo ''

temp_prefix="${tmpdir}"/temp_prefix
temp_outfile_0="${temp_prefix}".temp_outfile_0.txt
temp_outfile_1="${temp_prefix}".temp_outfile_1.txt
temp_outfile_2="${temp_prefix}".temp_outfile_2.txt
temp_outfile_3="${temp_prefix}".temp_outfile_3.txt
temp_outfile_4="${temp_prefix}".temp_outfile_4.txt

echo ''
echo 'Add sample id.'
echo 'awk -v sample='"$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmpout_genes_exact_filtered_splicing_tsv}" '>' "${temp_outfile_0}"
echo ''

awk -v sample="$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmpout_genes_exact_filtered_splicing_tsv}" > "${temp_outfile_0}"

report_gene_column_names='Gene.refGene,Gene.'"$wgEncodeGencodeBasic"',CLINVAR_GENEINFO_gene_name,spliceai_gene,spliceaiIndels_gene,FANTOM5_TSS_gene'

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

temp_outfile_incidentalome_1="${temp_prefix}".temp_outfile_incidentalome_1.txt
temp_outfile_incidentalome_2="${temp_prefix}".temp_outfile_incidentalome_2.txt

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

col_confident=0
col_tensorflow=0
echo 'col_confident_result=head -n 1' $temp_outfile_incidentalome_2 '| sed -e s/\t/\n/g | grep -n confident_splice_variant_scores'
col_confident_result=`head -n 1 $temp_outfile_incidentalome_2 | sed -e 's/\t/\n/g' | grep -n confident_splice_variant_scores`
echo 'col_tensorflow_result=head -n 1' $temp_outfile_incidentalome_2 '| sed -e s/\t/\n/g | grep -n needs_spliceai_tensorflow_scores'
col_tensorflow_result=`head -n 1 $temp_outfile_incidentalome_2 | sed -e 's/\t/\n/g' | grep -n needs_spliceai_tensorflow_scores`
echo ''
IFS=':' read -r -a array <<< "$col_confident_result"
col_confident="${array[0]}"
IFS=':' read -r -a array <<< "$col_tensorflow_result"
col_tensorflow="${array[0]}"
echo 'col_confident =' $col_confident
echo 'col_tensorflow =' $col_tensorflow

##### Extract the true positive splicing variants.
echo ''
echo 'Extract the true positive splicing variants.'
echo ''

echo 'awk -v col_confident='$col_confident 'BEGIN {FS="\t";OFS="\t"} {if ((NR==1)||($col_confident=="1")) {print $0}}' $temp_outfile_incidentalome_2 '>' $outfile
awk -v col_confident=$col_confident 'BEGIN {FS="\t";OFS="\t"} {if ((NR==1)||($col_confident=="1")) {print $0}}' $temp_outfile_incidentalome_2 > $outfile
echo ''

##### Extract the splicing variants for which we need spliceai-tensorflow to be run.
echo ''
echo 'Extract the splicing variants for which we need spliceai-tensorflow to be run.'
echo ''

echo 'awk -v col_tensorflow='$col_tensorflow 'BEGIN {FS="\t";OFS="\t"} {if ((NR==1)||($col_tensorflow=="1")) {print $0}}' $temp_outfile_incidentalome_2 '>' $outfile2
awk -v col_tensorflow=$col_tensorflow 'BEGIN {FS="\t";OFS="\t"} {if ((NR==1)||($col_tensorflow=="1")) {print $0}}' $temp_outfile_incidentalome_2 > $outfile2
echo ''

echo ''
echo 'outfile:' $outfile
echo 'outfile2:' $outfile2
echo ''

echo 'Finished!'
echo ''



