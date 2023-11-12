#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -l jobfs=100GB
#PBS -N cnvnator
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

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

#module load cmake/3.16.2
export ROOTSYS=$rootsys_for_cnvnator
export PATH=$PATH:$ROOTSYS/bin
export SHLIB_PATH=$ROOTSYS/lib
export LD_LIBRARY_PATH=$ROOTSYS/lib

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outfile_prefix="${outfile%.gz}"
outfile_prefix="${outfile_prefix%.vcf}"

rm -rf "${outfile_prefix}".cnvnator_output.root

# Extract read mapping
echo 'Extract read mapping'
echo $sw'/CNVnator_gadi/CNVnator/cnvnator -root' "${outfile_prefix}"'.cnvnator_output.root -tree' "${infile}" '-chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \'
echo '    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM'
$sw/CNVnator_gadi/CNVnator/cnvnator -root "${outfile_prefix}".cnvnator_output.root -tree "${infile}" -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
    chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM
echo ''

# Generate histogram
echo 'Generate histogram'
echo $cnvnator_executable '-root' "${outfile_prefix}"'.cnvnator_output.root -his 1000 -fasta' "${ref_fasta}" '-d' "${ref_fasta_chroms}"
$cnvnator_executable -root "${outfile_prefix}".cnvnator_output.root -his 1000 -fasta "${ref_fasta}" -d "${ref_fasta_chroms}"
echo ''

# Calculate statistics
echo 'Calculate statistics'
echo $cnvnator_executable '-root' "${outfile_prefix}"'.cnvnator_output.root -stat 1000'
$cnvnator_executable -root "${outfile_prefix}".cnvnator_output.root -stat 1000
echo ''

# Partition
echo 'Partition'
echo $cnvnator_executable '-root' "${outfile_prefix}"'.cnvnator_output.root -partition 1000'
$cnvnator_executable -root "${outfile_prefix}".cnvnator_output.root -partition 1000
echo ''

# Call CNVs
echo 'Call CNVs'
echo $cnvnator_executable '-root' "${outfile_prefix}"'.cnvnator_output.root -call 1000 >' "${outfile}"
$cnvnator_executable -root "${outfile_prefix}".cnvnator_output.root -call 1000 > "${outfile}"
echo ''

# https://github.com/abyzovlab/CNVnator/blob/master/README.md
# The output columns are as follows:
# CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0
# where,
# normalized_RD -- read depth normalized to 1.
# e-val1 -- is calculated using t-test statistics.
# e-val2 -- is from the probability of RD values within the region to be in the tails of a gaussian distribution describing frequencies of RD values in bins.
# e-val3 -- same as e-val1 but for the middle of CNV
# e-val4 -- same as e-val2 but for the middle of CNV
# q0 -- fraction of reads mapped with q0 quality

echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

