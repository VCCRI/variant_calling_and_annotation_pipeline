#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=30GB
#PBS -N bcftoolsExtractSample
#PBS -lstorage=gdata/abcd

set -euo pipefail

#infile=$1
#outdir=$2
#outfile=$3
#sw_and_refs=$4

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
module load bcftools

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

infile_basename=$(basename $infile)
tmpfile="${tmpdir}"/"${infile_basename}"
tmpfile_gz="${tmpdir}"/"${infile_basename}".gz

echo 'cp' $infile $tmpfile
cp $infile $tmpfile
echo ''

echo 'bgzip' $tmpfile
bgzip $tmpfile
echo ''

echo 'tabix -p vcf' $tmpfile_gz
tabix -p vcf $tmpfile_gz
echo ''

for sample in `bcftools query -l $infile`; do

  actual_outfile="${outfile/ExtractSamples/$sample}"
  tmpfile2="${actual_outfile%.vcf}".temp_file_needs_filtering.vcf

  #if [[ -f "$actual_outfile" ]]; then
  #  echo $actual_outfile 'already exists, no need to extract it again'
  #else

    echo 'bcftools view -s' $sample '-o' $tmpfile2 $tmpfile_gz
    bcftools view -s $sample -o $tmpfile2 $tmpfile_gz
    echo ''

    echo 'grep -P -v' '\t0\/0' $tmpfile2 '| grep -P -v' '\t\.\/\.' '| grep -P -v' '\t0\/\.' '| grep -P -v' '\t\.\/0' '>' $actual_outfile
    grep -P -v '\t0\/0' $tmpfile2 | grep -P -v '\t\.\/\.' | grep -P -v '\t0\/\.' | grep -P -v '\t\.\/0' > $actual_outfile
    echo ''
    #rm $tmpfile2

    echo 'Output file is:' $actual_outfile
    echo ''
  #fi

done

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"
