#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=30:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=2G
#PBS -N decomposePlatypus
#PBS -lstorage=gdata/abcd

set -euo pipefail

#infile=$1
#outdir=$2
#outfile=$3
#chrom=$4
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

module load htslib

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# Platypus headers have the following:
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">
##FORMAT=<ID=GOF,Number=.,Type=Float,Description="Goodness of fit value">
##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">
##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
#
# This means that NR and NV fields do not get their values split by vt decompose,
# so that multiple allelic values remain even when the alt has only one value (vt decompose converts multi-allelic to multiple single-allele rows)
# This would cause problems downstream in vpot when vpot needs to use NR and NV to calculate depth.
#
# Instead they should be the following:
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality as phred score">
##FORMAT=<ID=GOF,Number=1,Type=Float,Description="Goodness of fit value">
##FORMAT=<ID=NR,Number=A,Type=Integer,Description="Number of reads covering variant location in this sample">
##FORMAT=<ID=GL,Number=R,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">
##FORMAT=<ID=NV,Number=A,Type=Integer,Description="Number of reads containing variant in this sample">

infile_basename=$(basename $infile)
tmp_infile="${tmpdir}"/"${infile_basename}"

echo 'sed' 's/^##FORMAT=<ID=GQ,Number=.,/##FORMAT=<ID=GQ,Number=1,/' $infile '| sed' 's/^##FORMAT=<ID=GOF,Number=1,/##FORMAT=<ID=GOF,Number=.,/' '| sed' 's/^##FORMAT=<ID=NR,Number=.,/##FORMAT=<ID=NR,Number=A,/' '| sed' 's/^##FORMAT=<ID=GL,Number=.,/##FORMAT=<ID=GL,Number=G,/' '| sed' 's/^##FORMAT=<ID=NV,Number=.,/##FORMAT=<ID=NV,Number=A,/' '>' $tmp_infile
sed 's/^##FORMAT=<ID=GQ,Number=.,/##FORMAT=<ID=GQ,Number=1,/' $infile | sed 's/^##FORMAT=<ID=GOF,Number=1,/##FORMAT=<ID=GOF,Number=.,/' | sed 's/^##FORMAT=<ID=NR,Number=.,/##FORMAT=<ID=NR,Number=A,/' | sed 's/^##FORMAT=<ID=GL,Number=.,/##FORMAT=<ID=GL,Number=G,/' | sed 's/^##FORMAT=<ID=NV,Number=.,/##FORMAT=<ID=NV,Number=A,/' > $tmp_infile
echo ''

echo ''
echo 'Decompose' $tmp_infile 'to split multi-allelic variants to multiple single allelic variants, to produce output file' $outfile
echo ''

echo $sw'/vt/vt decompose -o' "$outfile" '-s' "$tmp_infile"
$sw/vt/vt decompose -o "$outfile" -s "$tmp_infile"
echo ''

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''
