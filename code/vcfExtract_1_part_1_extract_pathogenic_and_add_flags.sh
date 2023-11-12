#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=50G
#PBS -l ncpus=1
#PBS -l jobfs=100G
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

module load htslib
module load bcftools
module load python3
module load R/3.6.1
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
outfile_genes_patho_flags="${outdir}"/"${outfile_basename_prefix}"."${cohort}"_possiblyPathogenic.tsv

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

genes=""
regions=""
currdir=$(dirname $sw_and_refs)

if [[ "$cohort" == "CHD" ]]; then
  regions="${currdir}"/list_genes_CHD_955_genes_hg38_RefSeq_regions_colon_dash_format.txt
  genes="${currdir}"/list_genes_CHD_955_genes_withoutBars.txt
fi

if [[ "$cohort" == "DCM" ]]; then
  regions="${currdir}"/list_genes_DCM_all2022jan_909_genes_hg38_RefSeq_regions_colon_dash_format.txt
  genes="${currdir}"/list_genes_DCM_all2022jan_909_genes_withoutBars.txt
fi

# This is the final output file.
# This will continually have chromosome results added to it. 
# It will have the final variants and all the flags added.
rm -rf "${outfile_genes_patho_flags}"

##### List all the genes for a first rough pass.

echo ''
echo 'Extract variants from' $infile_prefix'*'$infile_suffix 'having input genes and input regions.'
echo 'Regions are' $regions 'Genes are' $genes
#echo 'Also extract the possibly_pathogenic variants from' $infile_prefix'*'$infile_suffix
echo ''

# Selecting by gene name can bring back too many hits, so select by gene name in context of tsv file.
# 
#	Gene.refGene=ABCB5				INFO_Gene.refGene
#	Gene.refGene=LINC02564\x3bLOC102723376
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

genes_string_for_filtering_vcf="${tmpdir}"/genes_string_for_filtering_vcf.txt
echo -e '^#' > $genes_string_for_filtering_vcf
#genes_string_for_filtering_tsv='^#'
while IFS= read -r inline; do

  echo -e "=${inline}" >> $genes_string_for_filtering_vcf
  echo -e "\\\x3b${inline}" >> $genes_string_for_filtering_vcf # annovar can produce \x3b in between gene name, which is encoding for ; This will pick up \x3b (correctly) and x3b (incorrectly)
  echo -e "|${inline}" >> $genes_string_for_filtering_vcf
  echo -e ";${inline}" >> $genes_string_for_filtering_vcf
  echo -e ":${inline}" >> $genes_string_for_filtering_vcf
  echo -e ",${inline}" >> $genes_string_for_filtering_vcf
  echo -e "|${inline}" >> $genes_string_for_filtering_vcf
  echo -e "_${inline}" >> $genes_string_for_filtering_vcf

  echo -e "${inline}\\t" >> $genes_string_for_filtering_vcf
  echo -e "${inline}\\\x3b" >> $genes_string_for_filtering_vcf # annovar can produce \x3b in between gene name, which is encoding for ; This will pick up \x3b (correctly) and x3b (incorrectly)
  echo -e "${inline};" >> $genes_string_for_filtering_vcf
  echo -e "${inline}:" >> $genes_string_for_filtering_vcf
  echo -e "${inline}," >> $genes_string_for_filtering_vcf
  echo -e "${inline}|" >> $genes_string_for_filtering_vcf

done < "$genes"

##### Roughly subset the input file to the genes and regions of interest.
##### Process by chromosomes because it gets too big for Rscript memory otherwise.

tmp_outfile_rough_genes_one_chrom="${tmpdir}"/tmp_outfile_rough_genes_one_chrom.tsv
tmp_outfile_genes_patho_one_chrom="${tmpdir}"/tmp_outfile_genes_patho_one_chrom.tsv
tmp_outfile_rough_genes="${tmpdir}"/tmp_outfile_rough_genes.tsv
tmp_outfile_genes_patho="${tmpdir}"/tmp_outfile_genes_patho.tsv
tmp_outfile_genes_patho_prev_chroms="${tmpdir}"/tmp_outfile_genes_patho_prev_chroms.tsv

mito="M"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
fi

#for chrom in "1" ; do
for chrom in {1..22} "$mito" "X" "Y" ; do

  echo ''
  echo 'Processing chrom' $chrom
  echo ''

  infile_chrom_vcf="${infile_prefix}${chrom}${infile_suffix}"
  infile_chrom_vcf_basename=$(basename $infile_chrom_vcf)
  infile_chrom_vcf_basename_prefix="${infile_chrom_vcf_basename%.gz}"
  infile_chrom_vcf_basename_prefix="${infile_chrom_vcf_basename_prefix%.vcf}"
  infile_chrom_vcf_prefix_basename=$(basename $infile_chrom_vcf_basename_prefix)
  infile_chrom_vcf_basename_prefix="${sample}".chrom"${chrom}"."${cohort}"

  tmpout_chrom_genes_vcf="${tmpdir}"/tmpout_chrom_genes_vcf.vcf
  tmpout_chrom_genes_tsv="${tmpdir}"/tmpout_chrom_genes_tsv.tsv
  tmpout_chrom_genes_refiltered_tsv="${tmpdir}"/tmpout_chrom_genes_refiltered_tsv.tsv
  tmpout_chrom_regions_vcf="${tmpdir}"/tmpout_chrom_regions_vcf.vcf
  tmpout_chrom_regions_tsv="${tmpdir}"/tmpout_chrom_regions_tsv.tsv
  tmpout_chrom_tier1tier2_tsv="${tmpdir}"/tmpout_chrom_tier1tier2_tsv.tsv
  tmpout_chrom_roughly_patho_vcf="${tmpdir}"/tmpout_chrom_roughly_patho_vcf.vcf
  tmpout_chrom_roughly_patho_tsv="${tmpdir}"/tmpout_chrom_roughly_patho_tsv.tsv

  echo ''
  echo 'Chrom:' $chrom 'Extract variants from' $infile_chrom_vcf 'having input genes, to produce' $tmpout_chrom_genes_vcf 'and' $tmpout_chrom_genes_tsv
  echo ''

  echo ''
  echo 'Chrom:' $chrom 'Extract variants by gene name.'
  echo 'grep -f' $genes_string_for_filtering_vcf $infile_chrom_vcf '>' $tmpout_chrom_genes_vcf
  echo ''

  set +e
  set +o pipefail
  grep -f "$genes_string_for_filtering_vcf" "${infile_chrom_vcf}" > "${tmpout_chrom_genes_vcf}"
  set -e
  set -e pipefail

  echo ''
  echo 'Chrom:' $chrom 'Extracted variants by gene name: convert vcf to tsv.'
  echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i' "${tmpout_chrom_genes_vcf}" '-o' "${tmpout_chrom_genes_tsv}" '-end_id YES'
  echo ''

  python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i "${tmpout_chrom_genes_vcf}" -o "${tmpout_chrom_genes_tsv}" -end_id YES

  # Don't refilter with only tabs on either side of gene name because there are gene fields that have more than one gene or value in the field, 
  # separated by various separators including underscore and \x3d and colon.
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
  #
  #echo ''
  #echo 'grep -P' $genes_string_for_filtering_tsv $tmpout_chrom_genes_tsv '>' $tmpout_chrom_genes_refiltered_tsv '|| True'
  #echo ''
  #
  #grep -P $genes_string_for_filtering_tsv $tmpout_chrom_genes_tsv > $tmpout_chrom_genes_refiltered_tsv || True

  tmp_infile_chrom_vcf="${tmpdir}"/"${infile_chrom_vcf_basename}"
  tmp_infile_chrom_vcf_gz="${tmpdir}"/"${infile_chrom_vcf_basename}".gz
  tmp_infile_chrom_vcf_gz_tbi="${tmpdir}"/"${infile_chrom_vcf_basename}".gz.tbi
  if [[ -f "${tmp_infile_chrom_vcf_gz_tbi}" ]]; then
    :
  else
    echo 'Chrom:' $chrom 'Copy and compress and index input vcf so that regions can be extracted by tabix.'
    echo 'Compress the input vcf' $infile_chrom_vcf 'to become' $tmp_infile_chrom_vcf_gz 'so that regions can be extracted by tabix'
    cp "${infile_chrom_vcf}" "${tmp_infile_chrom_vcf}"
    bgzip -f "${tmp_infile_chrom_vcf}"
    echo 'Index the copied and compressed input vcf' $tmp_infile_chrom_vcf_gz 'so that regions can be extracted by tabix'
    tabix -p vcf "${tmp_infile_chrom_vcf_gz}"
  fi

  echo ''
  echo 'Chrom:' $chrom 'Extract variants by region.'
  echo 'Use tabix to extract input_regions variants from' $tmp_infile_chrom_vcf_gz 'to produce' $tmpout_chrom_regions_vcf
  echo ''

  tabix -H "${tmp_infile_chrom_vcf_gz}" > "${tmpout_chrom_regions_vcf}"
  while IFS= read -r inline; do
    tabix "${tmp_infile_chrom_vcf_gz}" "${inline}" >> "${tmpout_chrom_regions_vcf}"
  done < "$regions"

  echo ''
  echo 'Chrom:' $chrom 'Extracted variants by region: convert vcf to tsv.'
  echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i' "${tmpout_chrom_regions_vcf}" '-o' "${tmpout_chrom_regions_tsv}" '-end_id YES'
  echo ''

  python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i "${tmpout_chrom_regions_vcf}" -o "${tmpout_chrom_regions_tsv}" -end_id YES

  echo ''
  echo 'Chrom:' $chrom 'Extracted variants by gene name and by region: concatenate these two tsv files to one file of unique tsv variants.'
  echo ''

  #echo 'head -n 1' "${tmpout_chrom_genes_tsv}" '>' "${tmp_outfile_rough_genes_one_chrom}"
  #head -n 1 "${tmpout_chrom_genes_tsv}" > "${tmp_outfile_rough_genes_one_chrom}"
  #echo 'cat' "${tmpout_chrom_genes_refiltered_tsv}" "${tmpout_chrom_regions_tsv}" '| grep -v "^CHROM" | sort | uniq >>' "${tmp_outfile_rough_genes_one_chrom}"
  #cat "${tmpout_chrom_genes_refiltered_tsv}" "${tmpout_chrom_regions_tsv}" | grep -v "^CHROM" | sort | uniq >> "${tmp_outfile_rough_genes_one_chrom}" || true

  echo 'head -n 1' "${tmpout_chrom_genes_tsv}" '>' "${tmp_outfile_rough_genes_one_chrom}"
  head -n 1 "${tmpout_chrom_genes_tsv}" > "${tmp_outfile_rough_genes_one_chrom}"
  echo 'cat' "${tmpout_chrom_genes_tsv}" "${tmpout_chrom_regions_tsv}" '| grep -v "^CHROM" | sort | uniq >>' "${tmp_outfile_rough_genes_one_chrom}"
  cat "${tmpout_chrom_genes_tsv}" "${tmpout_chrom_regions_tsv}" | grep -v "^CHROM" | sort | uniq >> "${tmp_outfile_rough_genes_one_chrom}" || true

  ##### Extract possibly pathogenic variants in the genes of interest, using the rough subset as input

  # Can't call this R script in one big hit. Too many records to match with too many genes. Have to break it up into subsets.

  echo ''
  echo 'Chrom:' $chrom 'Run the rough extract of variants by gene name and by region through extract_pathogenic_by_gene_name.R to get only possibly pathogenic variants in only matching genes.'
  echo 'Do the equivalent of the following command in subsets.'
  echo 'Rscript' "${sw}"'/victorchang_scripts/extract_pathogenic_by_gene_name.R' "${tmp_outfile_rough_genes_one_chrom}" "${tmp_outfile_genes_patho_one_chrom}" "${genes}"
  echo ''

  #Rscript "${sw}"/victorchang_scripts/extract_pathogenic_by_gene_name.R "${tmp_outfile_rough_genes_one_chrom}" "${tmp_outfile_genes_patho_one_chrom}" "${genes}"

  num_rows=$(wc -l $tmp_outfile_rough_genes_one_chrom)
  IFS=' ' read -r -a array <<< "$num_rows"
  num_rows="${array[0]}"
  num_loops=$(( num_rows / 5000 ))
  num_loops_extra=$(( num_rows % 5000 ))
  if [[ "$num_loops_extra" -gt 0 ]]; then
    num_loops=$(( num_loops + 1 ))
  fi

  :>$tmp_outfile_genes_patho_one_chrom
  for (( subset=1; subset<=$num_loops; subset++ )); do
    if [[ "$subset" -eq "$num_loops" ]]; then
      head_count=$num_rows
      tail_count=$num_loops_extra
    else
      head_count=$(( subset * 5000 ))
      tail_count=5000
    fi

    tmp_outfile_rough_genes_one_chrom_subset="${tmp_outfile_rough_genes_one_chrom}".subset.tsv
    tmp_outfile_rough_genes_one_chrom_hdr="${tmp_outfile_rough_genes_one_chrom}".hdr.tsv
    tmp_outfile_genes_patho_one_chrom_subset="${tmp_outfile_genes_patho_one_chrom}".subset.tsv
    tmp_outfile_genes_patho_one_chrom_hdr="${tmp_outfile_genes_patho_one_chrom}".hdr.tsv
    tmp_outfile_genes_patho_one_chrom_new="${tmp_outfile_genes_patho_one_chrom}".new.tsv

    echo ''
    echo 'Chrom:' $chrom 'Subset:' $subset 'Take a subset of input for feeding to extract_pathogenic_by_gene_name.R.'
    echo 'head -n 1' $tmp_outfile_rough_genes_one_chrom '>' $tmp_outfile_rough_genes_one_chrom_hdr
    head -n 1 $tmp_outfile_rough_genes_one_chrom > $tmp_outfile_rough_genes_one_chrom_hdr
    echo 'head -n' $head_count $tmp_outfile_rough_genes_one_chrom '| tail -n' $tail_count '| grep -v ^CHROM | cat' $tmp_outfile_rough_genes_one_chrom_hdr '- >' $tmp_outfile_rough_genes_one_chrom_subset
    head -n $head_count $tmp_outfile_rough_genes_one_chrom | tail -n $tail_count | grep -v '^CHROM' | cat $tmp_outfile_rough_genes_one_chrom_hdr - > $tmp_outfile_rough_genes_one_chrom_subset || true
    echo ''

    echo ''
    echo 'Chrom:' $chrom 'Subset:' $subset 'Run the rough extract of variants by gene name and by region through extract_pathogenic_by_gene_name.R to get only possibly pathogenic variants in only matching genes.'
    echo 'Rscript' "${sw}"'/victorchang_scripts/extract_pathogenic_by_gene_name.R' "${tmp_outfile_rough_genes_one_chrom_subset}" "${tmp_outfile_genes_patho_one_chrom_subset}" "${genes}"
    echo ''
    Rscript "${sw}"/victorchang_scripts/extract_pathogenic_by_gene_name.R "${tmp_outfile_rough_genes_one_chrom_subset}" "${tmp_outfile_genes_patho_one_chrom_subset}" "${genes}"

    # If the output is only a header, then don't use it. It may be the input header without any output header fields added to it.
    # For output, we want the true output header with all the output header fields.

    num_output_rows=$(wc -l $tmp_outfile_genes_patho_one_chrom_subset)
    IFS=' ' read -r -a array <<< "$num_output_rows"
    num_output_rows="${array[0]}"
    echo ''
    echo 'Chrom:' $chrom 'Subset:' $subset 'num_output_row' $num_output_rows 'in tmp_outfile_genes_patho_one_chrom_subset' $tmp_outfile_genes_patho_one_chrom_subset
    echo ''
    if [[ "$num_output_rows" -gt 1 ]]; then

      echo ''
      echo 'in if [[' "$num_output_rows" '-gt 1 ]]; then'
      echo 'Chrom:' $chrom 'Subset:' $subset 'Add subset of output from extract_pathogenic_by_gene_name.R to running output for this chromosome.'
      echo 'head -n 1' $tmp_outfile_genes_patho_one_chrom_subset '>' $tmp_outfile_genes_patho_one_chrom_hdr
      head -n 1 $tmp_outfile_genes_patho_one_chrom_subset > $tmp_outfile_genes_patho_one_chrom_hdr
      echo ''
      echo 'cat' $tmp_outfile_genes_patho_one_chrom $tmp_outfile_genes_patho_one_chrom_subset '| grep -v ^CHROM | cat' $tmp_outfile_genes_patho_one_chrom_hdr '- >' $tmp_outfile_genes_patho_one_chrom_new
      cat $tmp_outfile_genes_patho_one_chrom $tmp_outfile_genes_patho_one_chrom_subset | grep -v '^CHROM' | cat $tmp_outfile_genes_patho_one_chrom_hdr - > $tmp_outfile_genes_patho_one_chrom_new || true
      echo ''
      echo 'mv' $tmp_outfile_genes_patho_one_chrom_new $tmp_outfile_genes_patho_one_chrom
      mv $tmp_outfile_genes_patho_one_chrom_new $tmp_outfile_genes_patho_one_chrom
      echo ''
    fi

  done

  ##### Add sample id and gene flags and incidentalome flags

  infile_tsv_basename=$(basename $tmp_outfile_genes_patho_one_chrom)
  temp_prefix="${tmpdir}"/temp_prefix
  temp_outfile_0="${temp_prefix}".temp_outfile_0.txt
  temp_outfile_1="${temp_prefix}".temp_outfile_1.txt
  temp_outfile_2="${temp_prefix}".temp_outfile_2.txt
  temp_outfile_3="${temp_prefix}".temp_outfile_3.txt
  temp_outfile_4="${temp_prefix}".temp_outfile_4.txt

  echo ''
  echo 'Chrom:' $chrom 'Add sample id.'
  echo 'awk -v sample='"$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmp_outfile_genes_patho_one_chrom}" '>' "${temp_outfile_0}"
  echo ''

  awk -v sample="$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmp_outfile_genes_patho_one_chrom}" > "${temp_outfile_0}"

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

  if [[ -f "${outfile_genes_patho_flags}" ]]; then

    echo ''
    echo 'Chrom:' $chrom 'Make a copy of the possibly pathogenic variants in only matching genes of previous chromosomes.'
    echo 'cp' $outfile_genes_patho_flags $tmp_outfile_genes_patho_prev_chroms
    echo ''
    cp $outfile_genes_patho_flags $tmp_outfile_genes_patho_prev_chroms

    echo ''
    echo 'Chrom:' $chrom 'Concatenate the possibly pathogenic variants in only matching genes of previous chromosomes to the current chromosome just created.'
    echo 'grep -v "^CHROM"' $temp_outfile_incidentalome_2 '| cat' $tmp_outfile_genes_patho_prev_chroms '- >' $outfile_genes_patho_flags
    echo ''
    grep -v "^CHROM" $temp_outfile_incidentalome_2 | cat $tmp_outfile_genes_patho_prev_chroms - > $outfile_genes_patho_flags || true

  else

    echo ''
    echo 'Chrom:' $chrom 'This is the first chromosome so this results of possibly pathogenic variants in only matching genes is the currently running output.'
    echo 'cp' $temp_outfile_incidentalome_2 $outfile_genes_patho_flags
    echo ''
    cp $temp_outfile_incidentalome_2 $outfile_genes_patho_flags

  fi

done

echo ''
echo 'outfile_genes_patho_flags:' $outfile_genes_patho_flags
echo ''

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"



