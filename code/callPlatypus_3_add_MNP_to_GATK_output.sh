#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=30:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -N mergeGatkPlatypus
#PBS -lstorage=gdata/abcd

set -euo pipefail

#infile1=$1
#infile2=$2
#outdir=$3
#outfile=$4
#chrom=$5
#sw_and_refs=$6

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

module load htslib

tmp_hdr="${tmpdir}"/tmp_hdr.chrom"${chrom}".txt
tmp_hdr_platypus="${tmpdir}"/tmp_hdr_platypus.chrom"${chrom}".txt
tmp_hdr_platypus_samples="${tmpdir}"/tmp_hdr_platypus_samples.chrom"${chrom}".txt
tmp_hdr_gatk="${tmpdir}"/tmp_hdr_gatk.chrom"${chrom}".txt
tmp_hdr_gatk_samples="${tmpdir}"/tmp_hdr_gatk_samples.chrom"${chrom}".txt
tmp_platypus="${tmpdir}"/tmp_platypus.chrom"${chrom}".txt
tmp_platypus_with_mnvs="${tmpdir}"/tmp_platypus_with_mnvs.chrom"${chrom}".txt


tmp_platypus_uniq="${tmpdir}"/tmp_platypus_uniq.chrom"${chrom}".txt
tmp_gatk="${tmpdir}"/tmp_gatk.chrom"${chrom}".txt
tmp_gatk_cols1to5="${tmpdir}"/tmp_gatk_cols1to5.chrom"${chrom}".txt
tmp_gatk_rest="${tmpdir}"/tmp_gatk_rest.chrom"${chrom}".vcf
tmp_gatk_rest_minus_MNPs="${tmpdir}"/tmp_gatk_rest_minus_MNPs.chrom"${chrom}".vcf
tmp_body="${tmpdir}"/tmp_body.chrom"${chrom}".txt

grep '^##' $infile1 | grep -P 'INFO|FORMAT|FILTER' | grep -v '^##FORMAT=<ID=GT' | grep -v '^##FORMAT=<ID=GQ' > $tmp_hdr_platypus
echo -e '##INFO=<ID=IS_IT_MNP,Number=1,Type=String,Description="Does this location participate in a multinucleotide polymorphism. Values: IS_MNP">' >> $tmp_hdr_platypus
grep '^#CHROM' $infile1 > $tmp_hdr_platypus_samples

tabix -H $infile2 | grep -v '^#CHROM' > $tmp_hdr_gatk
tabix -H $infile2 | grep '^#CHROM' > $tmp_hdr_gatk_samples

list_platypus_samples=$(<$tmp_hdr_platypus_samples)
list_gatk_samples=$(<$tmp_hdr_gatk_samples)

if [[ $list_platypus_samples != $list_gatk_samples ]]; then
  echo 'Platypus and GATK samples are not the same'
  echo 'Platypus:'
  echo $list_platypus_samples
  echo 'GATK:'
  echo $list_gatk_samples
  echo 'EXITING'
  exit
fi

echo 'cat' $tmp_hdr_gatk $tmp_hdr_platypus $tmp_hdr_gatk_samples '>' $tmp_hdr
cat $tmp_hdr_gatk $tmp_hdr_platypus $tmp_hdr_gatk_samples > $tmp_hdr
echo ''

# Extract the Platypus MNVs - that is, variants that have 2 or more changed nucletides, not necessarily side by side
echo 'grep -v ^#' $infile1 '| awk BEGIN {FS="\t";OFS="\t"} {if ((length($4)>1) && (length($5)>1)) {printf $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 ";IS_IT_MNP=IS_MNP"; for (i=9; i<=NF; i++) {printf OFS $i}; printf "\n"}} >' $tmp_platypus
grep -v '^#' $infile1 | awk 'BEGIN {FS="\t";OFS="\t"} {if ((length($4)>1) && (length($5)>1)) {printf $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 ";IS_IT_MNP=IS_MNP"; for (i=9; i<=NF; i++) {printf OFS $i}; printf "\n"}}' > $tmp_platypus
echo ''

echo 'tabix' $infile2 'coords | awk BEGIN {FS="\t";OFS="\t"} {printf $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 ";IS_IT_MNP=IS_MNP"; for (i=9; i<=NF; i++) {printf OFS $i}; printf "\n"} >>' $tmp_gatk

:>$tmp_gatk
while IFS= read -r inline; do

  #echo $inline # this is the platypus MNP
  IFS=$'\t' read -r -a array <<< "$inline"
  chrom="${array[0]}"
  pos="${array[1]}"
  ref="${array[3]}"
  alt="${array[4]}"
  ref_size=${#ref}
  alt_size=${#alt}
  end=$(( pos + ref_size - 1 ))
  coords="${chrom}:${pos}-${end}"
  IFS=
  # these are the gatk SNPs that make up the platypus MNPs
  #echo 'tabix' $infile2 $coords '| awk BEGIN {FS="\t";OFS="\t"} {printf $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 ";IS_IT_MNP=IS_MNP"; for (i=9; i<=NF; i++) {printf OFS $i}; printf "\n"} >>' $tmp_gatk
  tabix $infile2 $coords | awk 'BEGIN {FS="\t";OFS="\t"} {printf $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 ";IS_IT_MNP=IS_MNP"; for (i=9; i<=NF; i++) {printf OFS $i}; printf "\n"}' >> $tmp_gatk

done < "$tmp_platypus"
echo ''

echo 'cut -d$''\t' '-f1-5' $tmp_gatk '>' $tmp_gatk_cols1to5
cut -d$'\t' -f1-5 $tmp_gatk > $tmp_gatk_cols1to5
echo ''

echo 'cp' $infile2 "${tmp_gatk_rest}"'.gz'
cp $infile2 "${tmp_gatk_rest}".gz
echo 'gunzip -f' "${tmp_gatk_rest}"'.gz'
gunzip -f "${tmp_gatk_rest}".gz
echo ''

# Remove the SNPs in GATK VCF that overlap MNPs because they will be added in again with the IS_MNP flag
echo 'grep -v -f' $tmp_gatk_cols1to5 $tmp_gatk_rest '>' $tmp_gatk_rest_minus_MNPs
grep -v -f $tmp_gatk_cols1to5 $tmp_gatk_rest > $tmp_gatk_rest_minus_MNPs
echo ''

# Platypus variants can be long and contain multiple MNPs.
# Extract the variants that match the gnomad MNV annovar annotation tables:
# 2bp MNVs, 3bp MNV, and SrS (dist=2), SrrS (dist=3), etc. up to dist=10, so variant will match gnomad SNV frequencies.
# This script includes outputting variants whose REF and ALT are only 2bp or 3bp long, and subset variants of longer variants.
echo 'awk -f' $sw'/victorchang_scripts/convert_long_variants_into_multiple_gnomad_MNV_format_variants.awk' $tmp_platypus '>' $tmp_platypus_with_mnvs
awk -f $sw/victorchang_scripts/convert_long_variants_into_multiple_gnomad_MNV_format_variants.awk $tmp_platypus > $tmp_platypus_with_mnvs
echo ''

# Get rid of any duplicates in our extracted smaller Platypus MNVs.
# We will not include the original longer Platypus MNVs because they cannot really be analysed because there will be no gnomad frequency annotations.
echo 'cat' $tmp_platypus_with_mnvs '| sort -k1,1V -k2,2n -k4,4 -k5,5 -k8,8 | uniq >' $tmp_platypus_uniq
cat $tmp_platypus_with_mnvs | sort -k1,1V -k2,2n -k4,4 -k5,5 -k8,8 | uniq > $tmp_platypus_uniq
echo ''

# Merge together the GATK vars, the SNPs in GATK VCF that overlap MNPs (have IS_MNP flag), and the Platypus MNPs (have IS_MNP flag)
echo 'cat' $tmp_gatk_rest_minus_MNPs $tmp_gatk $tmp_platypus_uniq '| grep -v ^# | sort -k1,1V -k2,2n -k4,4 -k5,5 -k8,8 | uniq >' $tmp_body
cat $tmp_gatk_rest_minus_MNPs $tmp_gatk $tmp_platypus_uniq | grep -v '^#' | sort -k1,1V -k2,2n -k4,4 -k5,5 -k8,8 | uniq > $tmp_body
echo ''

echo 'cat' $tmp_hdr $tmp_body '>' $outfile
cat $tmp_hdr $tmp_body > $outfile
echo ''

outfile2="${outfile%.gz}"
outfile2="${outfile2%.vcf}".removeAltAsterisk.vcf

echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if ($5 != "*") {print $0}}' $outfile '>' $outfile2
awk 'BEGIN {FS="\t";OFS="\t"} {if ($5 != "*") {print $0}}' $outfile > $outfile2
echo ''

echo 'outfile:' $outfile
echo 'outfile2:' $outfile2
echo ''
# bcftools view $outfile
# bcftools view -h $outfile

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''


