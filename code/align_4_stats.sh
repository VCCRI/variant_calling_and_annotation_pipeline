#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=8G
#PBS -l ncpus=1
#PBS -l jobfs=100G
#PBS -N stats
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

module load java/jdk-8.40
module load samtools

. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# Picard metrics

java -server -Xmx6g -Djava.io.tmpdir=$tmpdir \
    -jar $picard_jar CollectWgsMetrics \
    I="$infile" \
    O=$outfile \
    TMP_DIR=$tmpdir \
    R=$ref_fasta

# Samtools flagstat

samtools flagstat "${infile}" > "${outfile}".samtools_flagstat.txt

touch "${done_file}"
rm -f "${lock_file}"

echo 'Finished!'
echo ''
