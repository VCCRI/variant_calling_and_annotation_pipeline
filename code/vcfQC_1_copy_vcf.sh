#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l mem=4G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

# change first line of header in vcf file (4.2 becomes 4.0)
# GT 1/. causes plink --missing to crash, so change to 0/1
# Makes a copy of the input file
# Needed for certain software that doesn't deal with vcf 4.2 format

set -euo pipefail

#sample=$1
#infile=$2
#outdir=$3
#outfile=$4
#sw_and_refs=$5

# set environment variables for softwares and references
sw_and_refs=/my/directory/variant_calling_and_annotation_pipeline/code/where_are_softwares_and_references.sh
. "${sw_and_refs}"

sample="${batch}"
infile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".AllChrom.vcf.gz
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc
outfile1=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".AllChrom.GT01.vcf.gz
outfile2=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".AllChrom.GT01_VCFv4p0.vcf.gz
sw_and_refs=/my/directory/variant_calling_and_annotation_pipeline/code/where_are_softwares_and_references.sh

module load htslib
module load bcftools
module load samtools/1.10

. "${sw_and_refs}"

# A VCF that has been joint-called does not seem to cause plink to crash on genotypes 1/. and 0/. and thus no need to change those genotypes.
# A VCF created by merging singly-called samples can contain genotypes 1/. and 0/. that cause plink to crash and thus need to be changed to 0/1 and 0/1.
# Remove any ALT = * in case that is causing plink to crash

echo 'zcat' $infile '| awk to rmv ALT=* and rmv any Platypus MNV vars | bgzip >' $outfile1
{
    zcat $infile | awk 'BEGIN {FS="\t";OFS="\t"} {char1=substr($1,1,1); if ( (char1=="#") || (($5!="*") && ((length($4)==1) || (length($5)==1))) ) {print $0}}'
} | bgzip > $outfile1

echo 'tabix -p vcf' $outfile1
tabix -p vcf $outfile1
echo ''

echo 'zcat' $infile '| awk NR>1 | awk to rmv ALT=* | bgzip >' $outfile2
{
    echo "##fileformat=VCFv4.0"
    zcat $infile | awk 'NR>1' | awk 'BEGIN {FS="\t";OFS="\t"} {char1=substr($1,1,1); if ( (char1=="#") || (($5!="*") && ((length($4)==1) || (length($5)==1))) ) {print $0}}'
} | bgzip > $outfile2

echo 'tabix -p vcf' $outfile2
tabix -p vcf $outfile2
echo ''

echo 'outfile1' $outfile1
echo 'outfile2' $outfile2
echo ''
echo 'Finished!'
echo ''




