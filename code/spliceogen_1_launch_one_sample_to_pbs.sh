#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -l jobfs=20G
#PBS -N spliceogen
#PBS -lstorage=gdata/abcd

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

module load java/jdk-13.33
module load bedtools/

cd "${run_directory}"

./spliceogen_run_when_launched_in_pbs.sh -input "$input" -gtf "$gtf_for_spliceogen" -fasta "$fasta_for_spliceogen"

echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

