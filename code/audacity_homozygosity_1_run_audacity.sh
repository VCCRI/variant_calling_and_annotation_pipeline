#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N extracts
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile_prefix=$2
#infile_suffix=$3
#outdir=$4
#outfile=$5
#sw_and_refs=$6
#cohort=$7

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

module load intel-compiler/2021.1.1
module load R/4.0.0
module load htslib
module load bedtools
export PATH=$PATH:$sw/vcftools-0.14/bin/bin:/g/data/jb96/software/vcftools-0.14/bin:/g/data/jb96/software/vcftools-0.14
export PERL5LIB=$sw/vcftools-0.14/share/perl5:/g/data/jb96/software/vcftools-0.14/bin/share/perl5:/g/data/jb96/software/vcftools-0.14/src/perl:${PERL5LIB}

. "${sw_and_refs}"

outfile_basename=$(basename $outfile)

for chrom in {1..22}; do

  infile="${infile_prefix}${chrom}${infile_suffix}"

  sample_chrom_outdir="${outdir}"/"${sample}"_chr"${chrom}"
  mkdir -p $sample_chrom_outdir
  mkdir -p "${sample_chrom_outdir}"/"${sample}" # AUDACITY expects this sample-specific directory to be created already.

  out_label="${sample}"_chr"${chrom}"

  echo 'rm -rf' "${sample_chrom_outdir}"'/'"${out_label}"'.frq.gz'
  rm -rf "${sample_chrom_outdir}"/"${out_label}".frq.gz
  echo 'rm -rf' "${sample_chrom_outdir}"'/'"${out_label}"'.frq.gz.tbi'
  rm -rf "${sample_chrom_outdir}"/"${out_label}".frq.gz.tbi
  echo ''

  echo 'perl' $sw'/AUDACITY/AUDACITY/AUDACITYPrepare.pl -I' $infile '-F true -O' $sample_chrom_outdir '-L' $out_label
  perl $sw/AUDACITY/AUDACITY/AUDACITYPrepare.pl -I $infile -F true -O $sample_chrom_outdir -L $out_label
  echo ''

  echo 'perl' $sw'/AUDACITY/AUDACITY/AUDACITYAnalyze.pl -I' $sample_chrom_outdir '-L' $out_label '-F' "${out_label}"'.frq.gz -P1 0.1 -P2 0.1 -R2 0.001 -ND 100000 -NP 1 -AS' $genome_version '-O' $sample_chrom_outdir
  perl $sw/AUDACITY/AUDACITY/AUDACITYAnalyze.pl -I $sample_chrom_outdir -L $out_label -F "${out_label}".frq.gz -P1 0.1 -P2 0.1 -R2 0.001 -ND 100000 -NP 1 -AS $genome_version -O $sample_chrom_outdir
  echo ''

  echo 'ln -s' "${sample_chrom_outdir}"'/'"${sample}"'/'"${sample}"'_DIDOH3M2Regions.txt' "${sample_chrom_outdir}"'/'"${sample}"'/'"${sample}"'_DIDOH3M2Regions.txt.bed'
  ln -s "${sample_chrom_outdir}"/"${sample}"/"${sample}"_DIDOH3M2Regions.txt "${sample_chrom_outdir}"/"${sample}"/"${sample}"_DIDOH3M2Regions.txt.bed
  echo ''

done

echo ''
echo 'outfile:' "${sample_chrom_outdir}"/"${sample}"/"${sample}"_DIDOH3M2Regions.txt
echo ''

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

