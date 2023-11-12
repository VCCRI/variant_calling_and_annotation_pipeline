#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=2G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N splitVcfByChrom
#PBS -lstorage=scratch/abcd+gdata/abcd

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

exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "TERMINATING.  $1" >&2
        exit 1
    fi
}

module load bcftools
module load htslib

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo 'Split input' $infile 'by chromosome'
echo ''

infile_basename=$(basename $infile)
infile_basename_prefix="${infile_basename%.gz}"
infile_basename_prefix="${infile_basename_prefix%.vcf}"

outfile_part1="${outdir}"/"${infile_basename_prefix}".chrom
outfile_part2=.vcf

mito="M"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
fi

for chrom in {1..22} "$mito" "X" "Y" ; do

    out_vcf="${outfile_part1}${chrom}${outfile_part2}"

    echo ''
    echo 'Extract chromosome' $chrom 'of' $infile 'to produce chromosome output file' $out_vcf
    echo ''

    chr_chrom=chr"${chrom}"
    if [[ $genome_version == "hg19" ]]; then
      chr_chrom="${chrom}"
    fi

    echo 'bcftools view -r' ${chr_chrom} $infile '>' $out_vcf
    bcftools view -r ${chr_chrom} $infile > $out_vcf
    exit_on_error "Problem extracting vcf chrom: $chrom"

    echo ''
    echo 'Compress' $out_vcf
    echo ''

    bgzip -f "${out_vcf}"
    tabix -p vcf "${out_vcf}".gz
done

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

