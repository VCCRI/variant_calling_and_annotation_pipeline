#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=30:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N Platypus
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

exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "TERMINATING.  $1" >&2
        exit 1
    fi
}

module load python2-as-python
module load htslib

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

#tmpdir="${PBS_JOBFS}"/tmp
#mkdir -p "${tmpdir}"

echo 'python' $platypus 'callVariants \'
echo '  --bamFiles='$infile '--nCPU=3 \'
echo '  --refFile='$ref_fasta '\'
echo '  --output='$outfile '\'
echo '  --regions=$chrom --bufferSize=2000 --maxReads=25000000 --maxVariants=3'

python $platypus callVariants \
  --bamFiles=$infile --nCPU=3 \
  --refFile=$ref_fasta \
  --output=$outfile \
  --regions=$chrom --bufferSize=2000 --maxReads=25000000 --maxVariants=3

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''
