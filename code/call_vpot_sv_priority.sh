#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l mem=2G
#PBS -l ncpus=1
#PBS -l jobfs=2G
#PBS -N vpot
#PBS -lstorage=gdata/abcd

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

touch "${lock_file}"
rm -f "${queue_file}"

module load intel-mkl/2019.3.199
module load python3

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outfile_basename=$(basename $outfile)
outfile_basename_prefix="${outfile_basename%.tsv}"
outfile_basename_prefix="${outfile_basename_prefix%.txt}"
outfile_basename_prefix="${outfile_basename_prefix%.vcf}"
outfile_reformatted_to_call_vpot="${outdir}"/"${outfile_basename_prefix}"_reformatted_to_call_vpot.tsv

num_lines=$(wc -l $infile | cut -d' ' -f1)
if [[ $num_lines -eq 1 ]]; then
  echo ''
  echo 'There are no variants in' $infile
  echo ''
  echo 'Finished!'
  echo ''
  exit
fi

echo ''
echo 'sed -e s/\t\t/\t\.\t/g' "${infile}" '| sed -e s/\t\t/\t\.\t/g | sed -e s/\t$/\t\./g | awk -v sample='"${sample}" '-f' "${sw}"'/victorchang_scripts/reformat_sv_for_vpot.awk >' "${outfile_reformatted_to_call_vpot}"
echo ''
sed -e 's/\t\t/\t\.\t/g' "${infile}" | sed -e 's/\t\t/\t\.\t/g' | sed -e 's/\t$/\t\./g' | awk -v sample="${sample}" -f "${sw}"/victorchang_scripts/reformat_sv_for_vpot.awk > "${outfile_reformatted_to_call_vpot}"
echo ''

model=''
infile_list="${outdir}"/"${outfile_basename}".vpot_infile_list.txt
echo -e "${outfile_reformatted_to_call_vpot}\t${sample}" > "${infile_list}"

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

echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (($0 ~ /athogenic/) || ($0 ~ /p\.M1/) || ($2 != "0")) {print $0}}' "${vpot_output_file_renamed}" '| uniq >' "${vpot_output_file_vpot_gt_0}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (($0 ~ /athogenic/) || ($0 ~ /p\.M1/) || ($2 != "0")) {print $0}}' "${vpot_output_file_renamed}" | uniq > "${vpot_output_file_vpot_gt_0}"
echo ''

echo 'output file:' $vpot_output_file_vpot_gt_0
echo ''

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

