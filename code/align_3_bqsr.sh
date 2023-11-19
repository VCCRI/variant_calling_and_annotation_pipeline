#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=16G
#PBS -l ncpus=1
#PBS -l jobfs=16G
#PBS -N bqsr
#PBS -lstorage=gdata/abcd

set -euo pipefail

# BQSR --> Base Quality Score Recalibration
# Based on / adapted from GATK best practices
# NOTE < 1G jobfs used & ~ 10G mem ; longest job ~27hrs

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

. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# create BQSR table

out_bqsr_table="${outfile}".bqsr.table

$gatk_path --java-options "-server -Xms1g -Xmx14g -Djava.io.tmpdir=$tmpdir" BaseRecalibrator \
    -R "$ref_fasta" \
    -I "$infile" \
    --use-original-qualities \
    --known-sites "$dbsnp" \
    --known-sites "$known_indels" \
    --known-sites "$gold_std_indels" \
    -O "$out_bqsr_table"

# ApplyBQSR

$gatk_path --java-options "-server -Xms1g -Xmx14g -Djava.io.tmpdir=$tmpdir" ApplyBQSR \
    -R "$ref_fasta" \
    -I "$infile" \
    --bqsr-recal-file "$out_bqsr_table" \
    --static-quantized-quals 10 \
    --static-quantized-quals 20 \
    --static-quantized-quals 30 \
    --add-output-sam-program-record \
    --create-output-bam-index \
    --create-output-bam-md5 \
    --use-original-qualities \
    -O "$outfile"

# Validate the BAM file
java -server -Xms1g -Xmx14g -Djava.io.tmpdir=$tmpdir -jar $picard_jar \
      ValidateSamFile \
      I="$outfile" \
      TMP_DIR=$tmpdir \
      MODE=VERBOSE > "${outfile}".Picard_ValidateSamFile.txt

touch "${done_file}"
rm -f "${lock_file}"

echo 'Finished!'
echo ''
