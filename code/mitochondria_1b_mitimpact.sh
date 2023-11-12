#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N mitimpact
#PBS -lstorage=gdata/abcd

# Most sample chroms can process in walltime=5:00:00
# Some need more, eg. walltime=10:00:00

infile=$1
outdir=$2
outfile=$3
sw_and_refs=$4

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

module load python3/3.8.5
export PYTHONPATH=/my/directory/variant_calling_and_annotation_pipeline/victorchang_scripts/python3_3.8.5_packages_for_victorchang_scripts/lib/python3.8/site-packages

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo 'Call mitimpact - will need to call online https requests'
echo ''

echo 'python3' $sw'/victorchang_scripts/annotate_vcf_with_mitimpact_http_calls.py -i' $infile '-o' $outfile
python3 $sw/victorchang_scripts/annotate_vcf_with_mitimpact_http_calls.py -i $infile -o $outfile
echo ''

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'outfile' $outfile
echo ''
echo 'Finished!'
echo ''
