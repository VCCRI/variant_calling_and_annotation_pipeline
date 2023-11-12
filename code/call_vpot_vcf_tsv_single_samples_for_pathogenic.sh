#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l mem=2G
#PBS -l ncpus=1
#PBS -l jobfs=2G
#PBS -N vpot
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

set -euo pipefail

sample=$1
tool=$2
param_file=$3
infile=$4
outdir=$5
out_prefix=$6
sw_and_refs=$7

outfile=$out_prefix

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

#touch "${lock_file}"
#rm -f "${queue_file}"

module load intel-mkl/2019.3.199
module load python3

. "${sw_and_refs}"

rand_num="${sample}_call_vpot_vcf_tsv_single_samples_for_pathogenic"
tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outfile_basename=$(basename $outfile)
outfile_basename_prefix="${outfile_basename%.tsv}"
outfile_basename_prefix="${outfile_basename_prefix%.txt}"
outfile_basename_prefix="${outfile_basename_prefix%.vcf}"
outfile_reformatted_to_call_vpot="${outdir}"/"${outfile_basename_prefix}"_reformatted_to_call_vpot.tsv

echo ''
echo 'grep -v there_are_no_variants' "${infile}" '| sed -e s/\t\t/\t\.\t/g | sed -e s/\t\t/\t\.\t/g | sed -e s/\t$/\t\./g | awk -v sample='"${sample}" '-f' "${sw}"'/victorchang_scripts/reformat_vcf_tsv_for_vpot.awk >' "${outfile_reformatted_to_call_vpot}"
grep -v there_are_no_variants "${infile}" | sed -e 's/\t\t/\t\.\t/g' | sed -e 's/\t\t/\t\.\t/g' | sed -e 's/\t$/\t\./g' | awk -v sample="${sample}" -f "${sw}"/victorchang_scripts/reformat_vcf_tsv_for_vpot.awk > "${outfile_reformatted_to_call_vpot}"
echo ''

# Multinucleotide variants (MNPs) are annotated incorrectly and erroneously by annovar, so remove those annotations.
#
outfile_reformatted_to_call_vpot_rmvSomeAnnovar="${outdir}"/"${outfile_basename_prefix}"_reformatted_to_call_vpot.rmvSomeAnnovar.tsv
#
echo 'outfile_reformatted_to_call_vpot' $outfile_reformatted_to_call_vpot
c1=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n Func.refGene | cut -d':' -f1 | head -n 1)
c2=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n GeneDetail.refGene | cut -d':' -f1 | head -n 1)
c3=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n ExonicFunc.refGene | cut -d':' -f1 | head -n 1)
c4=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n AAChange.refGene | cut -d':' -f1 | head -n 1)
c5=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n Func.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c6=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n GeneDetail.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c7=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n ExonicFunc.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c8=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n AAChange.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c9=$(head -n 1 "${outfile_reformatted_to_call_vpot}" | sed -e 's/\t/\n/g' | grep -n VARIANT_TYPE | cut -d':' -f1 | head -n 1)
#
echo 'cat' $outfile_reformatted_to_call_vpot '| awk -v c1='"$c1" '-v c2='"$c2" '-v c3='"$c3" '-v c4='"$c4" '-v c5='"$c5" '-v c6='"$c6" '-v c7='"$c7" '-v c8='"$c8" '-v c9='"$c9" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0} else {if ((length($5)>1)&&(length($6)>1)) {$c1=".";$c2=".";$c3=".";$c4=".";$c5=".";$c6=".";$c6=".";$c7=".";$c8=".";$c9="."; print $0} else {print $0}}}' '>' $outfile_reformatted_to_call_vpot_rmvSomeAnnovar
cat $outfile_reformatted_to_call_vpot | awk -v c1="$c1" -v c2="$c2" -v c3="$c3" -v c4="$c4" -v c5="$c5" -v c6="$c6" -v c7="$c7" -v c8="$c8" -v c9="$c9" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0} else {if ((length($5)>1)&&(length($6)>1)) {$c1=".";$c2=".";$c3=".";$c4=".";$c5=".";$c6=".";$c7=".";$c8=".";$c9="."; print $0} else {print $0}}}' > $outfile_reformatted_to_call_vpot_rmvSomeAnnovar
echo ''

model=''
infile_list="${outdir}"/"${outfile_basename}".vpot_infile_list.txt
echo -e "${outfile_reformatted_to_call_vpot_rmvSomeAnnovar}\t${sample}" > "${infile_list}"
#echo -e "${outfile_reformatted_to_call_vpot}\t${sample}" > "${infile_list}"

echo "Input parameters to VPOT:"
echo "Tool option : "$tool
echo "Output dir + prefix : "$out_prefix
echo "input file : "$infile_list
echo "contents of input file : "$infile
echo "processing parameter : "$param_file
echo ''

echo ''
echo 'python3 '"${sw}"'/VPOT/VPOT.py' $tool $out_prefix $infile_list $param_file
echo ''
python3 "${sw}"/VPOT/VPOT.py $tool $out_prefix $infile_list $param_file
echo ''

vpot_output_file_pattern="${out_prefix}"final_output_file_*.txt
vpot_output_file_renamed="${out_prefix}"_final_output_file.txt
vpot_output_file_vpot_gt_0="${out_prefix}"_final_output_file.vpot_gt_0.txt

echo 'Rename the vpot output file to remove the random number in the file name.'
vpot_output_file=$(ls -t $vpot_output_file_pattern | head -n 1)
echo 'mv' $vpot_output_file $vpot_output_file_renamed
mv $vpot_output_file $vpot_output_file_renamed
echo ''

# Sort by chrom+pos,
# assign the same highest ranking and score to overlapping variants such as MNVs and their singleton SNPs,
# then go back to being sorted by score and ranking.
# Thus overlapping GATK/Platypus/long_MNPs/3bp_MNPs will have the same highest vpot score and thus appear as a cluster in the vpot report.

vpot_output_file_adjusted_scores="${tmpdir}"/vpot_output_file_adjusted_scores.txt
tmphdr="${tmpdir}"/tmphdr.txt
head -n 1 $vpot_output_file_renamed > $tmphdr

echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $vpot_output_file_renamed '| sort -k3,3V -k4,4V -k6,6V -k7,7V | \'
echo '  awk -f' $sw'/victorchang_scripts/adjust_vpot_score_for_MNVs.awk | sort -k2,2Vr -k1,1Vr -k3,3V -k4,4V -k6,6V -k7,7V | cat' $tmphdr '- >' $vpot_output_file_adjusted_scores
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $vpot_output_file_renamed | sort -k3,3V -k4,4V -k6,6V -k7,7V | \
  awk -f $sw/victorchang_scripts/adjust_vpot_score_for_MNVs.awk | sort -k2,2Vr -k1,1Vr -k3,3V -k4,4V -k6,6V -k7,7V | cat $tmphdr - > $vpot_output_file_adjusted_scores
echo ''

# Filter out variants having vpot ranking of zero.
# However, don't filter out variants that have pathogenic annotation (unfortunately will also catch conflicting...pathogenicity).

vep_impact_col=-1
cols_string=$(head -n 1 $vpot_output_file_adjusted_scores)
IFS=$'\t' read -r -a array <<< "$cols_string"
i=0
for this_col in "${array[@]}"; do
  i=$(( i + 1 ))
  if [[ $this_col == "VEP_IMPACT" ]]; then
    vep_impact_col=$i
  fi
  if [[ $this_col == "INFO_VEP_IMPACT" ]]; then
    vep_impact_col=$i
  fi
done

if [[ $vep_impact_col -gt -1 ]]; then
  echo 'awk -v vep_impact_col='"$vep_impact_col" 'BEGIN {FS="\t";OFS="\t"} {if (('$vep_impact_col' == "HIGH") || ($0 ~ /athogenic/) || ($0 ~ /p\.M1/) || ($2 != "0")) {print $0}}' "${vpot_output_file_adjusted_scores}" '| uniq >' "${vpot_output_file_vpot_gt_0}"
  awk -v vep_impact_col="$vep_impact_col" 'BEGIN {FS="\t";OFS="\t"} {if (($vep_impact_col == "HIGH") || ($0 ~ /athogenic/) || ($0 ~ /p\.M1/) || ($2 != "0")) {print $0}}' "${vpot_output_file_adjusted_scores}" | uniq > "${vpot_output_file_vpot_gt_0}"
else
  echo 'awk -v' 'BEGIN {FS="\t";OFS="\t"} {if (($0 ~ /athogenic/) || ($0 ~ /p\.M1/) || ($2 != "0")) {print $0}}' "${vpot_output_file_adjusted_scores}" '| uniq >' "${vpot_output_file_vpot_gt_0}"
  awk 'BEGIN {FS="\t";OFS="\t"} {if (($0 ~ /athogenic/) || ($0 ~ /p\.M1/) || ($2 != "0")) {print $0}}' "${vpot_output_file_adjusted_scores}" | uniq > "${vpot_output_file_vpot_gt_0}"
fi
echo ''

# Multinucleotide variants (MNPs) are annotated incorrectly and erroneously by annovar, so remove those annotations.
#
vpot_output_file_vpot_gt_0_rmvSomeAnnovar="${vpot_output_file_vpot_gt_0%.tsv}".rmvSomeAnnovar.tsv
#
c1=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n Func.refGene | cut -d':' -f1 | head -n 1)
c2=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n GeneDetail.refGene | cut -d':' -f1 | head -n 1)
c3=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n ExonicFunc.refGene | cut -d':' -f1 | head -n 1)
c4=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n AAChange.refGene | cut -d':' -f1 | head -n 1)
c5=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n Func.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c6=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n GeneDetail.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c7=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n ExonicFunc.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
c8=$(head -n 1 "${vpot_output_file_vpot_gt_0}" | sed -e 's/\t/\n/g' | grep -n AAChange.wgEncodeGencodeBasic | cut -d':' -f1 | head -n 1)
#
echo 'cat' $vpot_output_file_vpot_gt_0 '| awk -v c1='"$c1" '-v c2='"$c2" '-v c3='"$c3" '-v c4='"$c4" '-v c5='"$c5" '-v c6='"$c6" '-v c7='"$c7" '-v c8='"$c8" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0} else {if ((length($7)>1)&&(length($8)>1)) {$c1=".";$c2=".";$c3=".";$c4=".";$c5=".";$c6=".";$c6=".";$c7=".";$c8="."; print $0} else {print $0}}}' '>' $vpot_output_file_vpot_gt_0_rmvSomeAnnovar
cat $vpot_output_file_vpot_gt_0 | awk -v c1="$c1" -v c2="$c2" -v c3="$c3" -v c4="$c4" -v c5="$c5" -v c6="$c6" -v c7="$c7" -v c8="$c8" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0} else {if ((length($7)>1)&&(length($8)>1)) {$c1=".";$c2=".";$c3=".";$c4=".";$c5=".";$c6=".";$c7=".";$c8="."; print $0} else {print $0}}}' > $vpot_output_file_vpot_gt_0_rmvSomeAnnovar
echo ''

echo 'output file:' $vpot_output_file_vpot_gt_0
echo 'output file2:' $vpot_output_file_vpot_gt_0_rmvSomeAnnovar
echo ''

echo ''
echo 'Finished!'
echo ''; echo ''

touch "${done_file}"
rm -f "${lock_file}"

