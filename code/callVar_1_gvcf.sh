#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=6G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N HaplotypeCaller
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

exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "TERMINATING.  $1" >&2
        exit 1
    fi
}

module load python3-as-python
module load java/jdk-8.40

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo 'Run gatk HaplotypeCaller on' $infile 'to produce the output chromosome files that will later on be merged to produce' $outfile
echo ''

for chrom in $chroms ; do

    chr_chrom="chr${chrom}"
    if [[ $genome_version == "hg19" ]]; then
      chr_chrom="${chrom}"
    fi

    out_gvcf="${outfile_part1}${chrom}${outfile_part2}"

    echo ''
    echo 'Run gatk HaplotypeCaller on chromosome' $chrom 'of' $infile 'to produce chromosome output file' $out_gvcf
    echo ''

    $gatk_path --java-options "-server -Xms1g -Xmx4g -Djava.io.tmpdir=$tmpdir" \
        HaplotypeCaller \
        -R "$ref_fasta" \
        --dbsnp "$gatk_dbsnp" \
        -I "$infile" \
        -ERC GVCF \
        -O "$out_gvcf" \
        -L "${chr_chrom}"
    exit_on_error "Problem creating gvcf for chrom: $chrom"
done

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''
