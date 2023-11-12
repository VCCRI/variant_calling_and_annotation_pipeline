#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=6G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N compare
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile=$2
#outfile=$3
#database_file=$4
#sw_and_refs=$5
#include_low_qual_flag=$6
#exclude_samples=$7

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


# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo 'Rscript' $sw'/victorchang_scripts/compare_strucvars_with_other_samples.R' "${sample}" "${infile}" "${outfile}" "${database_file}" "${include_low_qual_flag}" "${exclude_samples}"

Rscript $sw/victorchang_scripts/compare_strucvars_with_other_samples.R "${sample}" "${infile}" "${outfile}" "${database_file}" "${include_low_qual_flag}" "${exclude_samples}"

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'

touch "${done_file}"
rm -f "${lock_file}"

