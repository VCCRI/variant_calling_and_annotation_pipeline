#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=64G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -N spliceaitensonflow
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile=$2
#infile2=$3
#outdir=$4
#outfile_basename=$5
#sw_and_refs=$6

# Set the job to terminate when fail, so that we don't get the Finished! message and we know it failed.
# However, the counts can fail, so don't terminate when they fail.
set -euo pipefail

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outfile="${outdir}"/"${outfile_basename}"
outfile_vcf="${outfile}".vcf
outfile_vcf_split_info="${outfile}".split_info.vcf
outfile_tsv="${outfile}".tsv

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

module load R/3.6.1
module unload intel-fc intel-cc
module load intel-compiler/2019.3.199

module load bcftools
module load cuda/10.1
module load python3
module load hdf5/1.10.5
export PYTHONPATH=$pythonpath_for_spliceai
module load tensorflow/2.0.0

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

echo ''
echo '### 1) Create a vcf file from the input tsv file'
echo ''

tmp_vcf="${tmpdir}"/"${sample}".tmp_vcf_for_tensorflow.vcf

echo 'grep -P' '^##file|^##FILTER=|^##contig=' $infile2 '>' $tmp_vcf
grep -P '^##file|^##FILTER=|^##contig=' $infile2 > $tmp_vcf
echo ''
echo '##INFO=<ID=SAMPLE_FOR_SPLICEAI,Number=1,Type=String,Description="Sample ID">' >> $tmp_vcf
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> $tmp_vcf
echo 'awk -v sample='"$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1,$2,$4,$5,$6,$7,$8, "SAMPLE_FOR_SPLICEAI=" sample }}' $infile '>>' $tmp_vcf
awk -v sample="$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1,$2,$4,$5,$6,$7,$8, "SAMPLE_FOR_SPLICEAI=" sample }}' $infile >> $tmp_vcf
echo ''

echo ''
echo '### 2) Submit vcf file to spliceai program (look 500 bp either side)'
echo ''

echo '"${sw}"/spliceai_pip_install/spliceai/bin/spliceai -D 500 -I' $tmp_vcf '-O' $outfile_vcf '-R' $spliceai_reference '-A grch38'
"${sw}"/spliceai_pip_install/spliceai/bin/spliceai -D 500 -I $tmp_vcf -O $outfile_vcf -R $spliceai_reference -A grch38
echo ''

count=`grep -v '^#' $outfile_vcf | grep -v '^CHROM' | grep -v '^chrom' | wc -l | cut -d' ' -f1 || true`
echo $count 'data records in' $outfile_vcf
echo ''

echo ''
echo '### 3) Split spliceai tensorflow vcf result into separate info fields in a vcf'
echo ''

echo 'python3' $sw'/victorchang_scripts/split_vcf_info_field_into_multiple_fields.py -i' $outfile_vcf '-o' $outfile_vcf_split_info '-f SpliceAI -d |'
python3 $sw/victorchang_scripts/split_vcf_info_field_into_multiple_fields.py -i $outfile_vcf -o $outfile_vcf_split_info -f SpliceAI -d '|'
echo ''

echo ''
echo '### 4) Convert spliceai tensorflow vcf result to tab-delimited'
echo ''

echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i' $outfile_vcf_split_info '-o' $outfile_tsv '-end_id YES'
python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i $outfile_vcf_split_info -o $outfile_tsv -end_id YES
echo ''

echo ''
echo 'output_vcf:' $outfile_vcf
echo 'output_tsv:' $outfile_tsv
echo ''

echo 'Finished!'


touch "${done_file}"
rm -f "${lock_file}"

