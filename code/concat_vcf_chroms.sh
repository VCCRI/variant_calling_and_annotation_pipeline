#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=7G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N concat_vcf_chroms
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile_prefix=$2
#infile_suffix=$3
#outdir=$4
#outfile=$5
#sw_and_refs=$6

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

touch "${lock_file}"
rm -f "${queue_file}"

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

infile_for_header="${infile_prefix}1${infile_suffix}"
infile_for_header_basename=$(basename $infile_for_header)
infile_header="${infile_for_header}".header.vcf
grep '^#' "${infile_for_header}" > "${infile_header}"

infile_dirname=$(dirname $infile_for_header)
infile_prefix_basename=$(basename $infile_prefix)
cd "${infile_dirname}"

mito="M"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
fi

cat_string=''
for chrom in {1..22} "$mito" "X" "Y" ; do

  echo 'Processing chrom' $chrom

  infile="${infile_prefix_basename}${chrom}${infile_suffix}"

  if [ "$cat_string" = "" ]; then
    cat_string="${infile}"
  else
    cat_string="${cat_string} ${infile}"
  fi

done

cat $cat_string | grep -v '^#' | cat "${infile_header}" - > "${outfile}"

rm "${infile_header}"

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished! Output file is' $outfile
echo ''


