#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=20G
#PBS -lstorage=gdata/abcd

# change first line of header in vcf file (4.2 becomes 4.0)
# Makes a copy of the input file
# Needed for certain software that doesn't deal with vcf 4.2 format

set -euo pipefail

# set environment variables for softwares and references
sw_and_refs=/my/directory/variant_calling_and_annotation_pipeline/code/where_are_softwares_and_references.sh
. "${sw_and_refs}"

# Use the GATK VCF before normalisation, which is what is used for CHD BATCH_7
indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_vqsr
infile="${indir}"/"${batch}"__allChrom.vqsr.vcf.gz

outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc
outfile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".AllChrom.vcf
outfile_gz="${outfile}".gz

module load htslib
module load bcftools
module load samtools

outfile_basename=$(basename $outfile)
tmpdir="${PBS_JOBFS}"/tmp/"${outfile_basename}"
mkdir -p "${tmpdir}"
tmpdir_for_bcftools="${PBS_JOBFS}"/tmpdir_for_bcftools
mkdir -p "${tmpdir_for_bcftools}"

echo 'bcftools sort' $infile '-O z -o' $outfile_gz '-T' $tmpdir_for_bcftools
bcftools sort $infile -O z -o $outfile_gz -T $tmpdir_for_bcftools
echo ''

echo 'tabix -p vcf' $outfile_gz
tabix -p vcf $outfile_gz
echo ''

echo 'Output file is' $outfile_gz
echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"



