#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=2G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N decompose
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
module load samtools
module load htslib

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

infile_basename=$(basename $infile)
infile_basename_prefix="${infile_basename%.gz}"
infile_basename_prefix="${infile_basename_prefix%.vcf}"
infile_basename_prefix="${infile_basename_prefix%.g}"

# After benchmarking this pipeline, change this job so as to not keep all these intermediate files, by creating them on PBS_JOBFS temporary file space.

decompose_vcf="${outdir}"/"${infile_basename_prefix}".decompose.vcf.gz
blocksub_vcf="${outdir}"/"${infile_basename_prefix}".decompose_blocksub.vcf.gz
normalize_vcf="${outdir}"/"${infile_basename_prefix}".decompose_blocksub_normalize.vcf.gz

# decompose
echo ''
echo 'Decompose' $infile 'to produce output file' $decompose_vcf
echo ''

$sw/vt/vt decompose -o "$decompose_vcf" -s "$infile"
tabix -p vcf "$decompose_vcf"

# blocksub
echo ''
echo 'Blocksub' $decompose_vcf 'to produce output file' $blocksub_vcf
echo ''

$sw/vt/vt decompose_blocksub -o "$blocksub_vcf" -a "$decompose_vcf"
tabix -p vcf "$blocksub_vcf"

# normalize
echo ''
echo 'Normalize' $blocksub_vcf 'to produce output file' $normalize_vcf
echo ''

$sw/vt/vt normalize -r $ref_fasta -o "$normalize_vcf" "$blocksub_vcf"
tabix -p vcf "$normalize_vcf"

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

