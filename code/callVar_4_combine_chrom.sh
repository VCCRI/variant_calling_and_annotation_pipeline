#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=01:00:00
#PBS -l mem=4G
#PBS -l ncpus=1
#PBS -l jobfs=16G
#PBS -N GatherVcfs
#PBS -lstorage=gdata/abcd

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

module load python3-as-python
module load java/jdk-8.40
module load htslib

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo 'Run gatk GatherVcfs' $infile_prefix $infile_suffix 'for chromosomes' $chroms 'to produce output file' $outfile
echo ''

IFS=':' read -r -a array <<< "$chroms"
input_string=""
for chrom in "${array[@]}"; do
    if [[ "$input_string" == "" ]]; then
        input_string="-I "$infile_prefix$chrom$infile_suffix
    else
        temp_string="${input_string} -I "$infile_prefix$chrom$infile_suffix
        input_string=$temp_string
    fi
done

${gatk_path} --java-options "-server -Xms1g -Xmx2g -Djava.io.tmpdir=$tmpdir" \
    GatherVcfs \
    $input_string \
    -O ${outfile}

echo ''
echo 'Create index for output VCF by running tabix -p vcf' $outfile
#echo 'Run gatk IndexFeatureFile for' $outfile
echo ''

tabix -p vcf "${outfile}"

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

