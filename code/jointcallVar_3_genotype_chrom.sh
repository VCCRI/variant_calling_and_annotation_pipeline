#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -N GenotypeGVCFs
#PBS -lstorage=gdata/abcd

# 10 hours walltime is not enough for chrom 1-7,X for trios

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

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo 'Run gatk GenotypeGVCFs on chromosome' $chrom 'of' $infile 'to produce output file' $outfile
echo ''

echo ${gatk_path} '--java-options "-server -Xms1g -Xmx5g -Djava.io.tmpdir=$tmpdir" \'
echo '    GenotypeGVCFs \'
echo '    -R' ${ref_fasta} '\'
echo '    -O' ${outfile} '\'
echo '    -D' ${gatk_dbsnp} '\'
echo '    -G StandardAnnotation \'
echo '    --only-output-calls-starting-in-intervals \'
echo '    --use-new-qual-calculator \'
echo '    --genomicsdb-shared-posixfs-optimizations \'
echo '    -V gendb://'${infile} '\'
echo '    -L' ${chrom}

${gatk_path} --java-options "-server -Xms1g -Xmx5g -Djava.io.tmpdir=$tmpdir" \
    GenotypeGVCFs \
    -R ${ref_fasta} \
    -O ${outfile} \
    -D ${gatk_dbsnp} \
    -G StandardAnnotation \
    --only-output-calls-starting-in-intervals \
    --use-new-qual-calculator \
    --genomicsdb-shared-posixfs-optimizations \
    -V gendb://${infile} \
    -L ${chrom}

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

