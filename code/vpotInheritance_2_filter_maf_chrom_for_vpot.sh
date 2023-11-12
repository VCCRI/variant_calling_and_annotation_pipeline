#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=60G
#PBS -l ncpus=1
#PBS -l jobfs=2G
#PBS -N vpot
#PBS -lstorage=gdata/abcd

#param_file=$1
#infile=$2
#outdir=$3
#outfile=$4
#sw_and_refs=$5

set -euo pipefail

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

module load intel-mkl/2019.3.199
module load python3

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

out_prefix="${outfile%.txt}"
outfile_basename=$(basename $outfile)

tmp_reformatted="${tmpdir}"/"${outfile_basename}".tmp_reformatted.vcf

echo 'awk -v min_vaf_for_snvs=0.1 -v min_depth_for_snvs=1 -v min_vaf_for_mnvs=0.2 -v min_depth_for_mnvs=10 -f' $sw'/victorchang_scripts/filter_vcf_snv_and_mnv_vars_for_min_depth_and_min_vaf.awk' $infile '| awk -f' $sw'/victorchang_scripts/reformat_vcf_for_vpot_and_set_NR_NV_to_max_values.awk | awk' 'BEGIN {FS="\t";OFS="\t"} {if ($5 != "*") {print $0}}' '>' $tmp_reformatted
awk -v min_vaf_for_snvs=0.1 -v min_depth_for_snvs=1 -v min_vaf_for_mnvs=0.2 -v min_depth_for_mnvs=10 -f $sw/victorchang_scripts/filter_vcf_snv_and_mnv_vars_for_min_depth_and_min_vaf.awk $infile | awk -f $sw/victorchang_scripts/reformat_vcf_for_vpot_and_set_NR_NV_to_max_values.awk | awk 'BEGIN {FS="\t";OFS="\t"} {if ($5 != "*") {print $0}}' > $tmp_reformatted
echo ''

list_of_samples=$(grep '^#CHROM' $infile)
echo 'list_of_samples' $list_of_samples
IFS=$'\t' read -r -a array <<< "$list_of_samples"

infile_list="${out_prefix}".vpot_infile_list.txt

:>$infile_list
num_cols="${#array[@]}"
for i in $(seq 10 1 $num_cols); do
  j=$(( i - 1 ))
  sample=${array[j]}
  echo -e "${tmp_reformatted}\t${sample}" >> "${infile_list}"
done

echo 'python3' $sw'/VPOT/VPOT.py priority' $out_prefix $infile_list $param_file
python3 $sw/VPOT/VPOT.py priority $out_prefix $infile_list $param_file
echo ''

vpot_output_file_pattern="${out_prefix}"final_output_file_*.txt
vpot_output_file_renamed="${out_prefix}"_final_output_file.txt

echo 'Rename the vpot output file to remove the random number in the file name.'
vpot_output_file=$(ls -t $vpot_output_file_pattern | head -n 1)
echo 'mv' $vpot_output_file $vpot_output_file_renamed
mv $vpot_output_file $vpot_output_file_renamed
echo ''

echo 'output:' $vpot_output_file_renamed
echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

