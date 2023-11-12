#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=8G
#PBS -l ncpus=1
#PBS -l jobfs=80G
#PBS -N stats
#PBS -lstorage=gdata/abcd

#sample=$1
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

module load java/jdk-8.40
module load samtools
module load R/3.6.1

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo 'Run picard.jar CollectInsertSizeMetrics' $picard_jar

java -server -Xmx6g -Djava.io.tmpdir=$tmpdir \
      -jar $picard_jar CollectInsertSizeMetrics \
      I="${infile}" \
      O="${outdir}"/"${sample}".picard_insert_size_metrics.txt \
      H="${outdir}"/"${sample}".picard_insert_size_histogram.pdf \
      M=0.5

echo ''
echo 'Run picard.jar CollectJumpingLibraryMetrics' $picard_jar

java -server -Xmx6g -Djava.io.tmpdir=$tmpdir \
      -jar $picard_jar CollectJumpingLibraryMetrics \
      I="${infile}" \
      O="${outdir}"/"${sample}".picard_jumping_metrics.txt

echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

