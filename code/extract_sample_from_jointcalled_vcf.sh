#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=10B
#PBS -N bcftoolsExtractSample
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

set -euo pipefail

#sample=$1
#infile=$2
#outdir=$3
#outfile=$4
#sw_and_refs=$5

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

module load htslib/1.9
module load bcftools/1.9

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo 'bcftools view -c1 -Oz -s' $sample '-o' $outfile $infile
bcftools view -c1 -Oz -s $sample -o $outfile $infile
echo ''

echo 'tabix -p vcf' $outfile
tabix -p vcf $outfile
echo ''

echo 'Finished! Output file is:' $outfile
echo ''

touch "${done_file}"
rm -f "${lock_file}"
