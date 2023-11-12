#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N concat
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

module load htslib/1.9

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

infile_prefix_basename=$(basename $infile_prefix)
infile_header="${infile_prefix_basename}"1.header.vcf

infile_for_header="${infile_prefix}"1"${infile_suffix}"
if [[ $infile_for_header == *.gz ]]; then
  infile_for_header_gz="${tmpdir}"/"${infile_prefix_basename}"1"${infile_suffix}"
  infile_for_header="${infile_for_header_gz%.gz}"
fi

mito="M"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
fi

cat_string=''
for chrom in {1..22} "$mito" "X" "Y" ; do

  echo 'chrom' $chrom

  temp_infile="${infile_prefix}${chrom}${infile_suffix}"
  temp_outfile="${tmpdir}"/"${infile_prefix_basename}${chrom}${infile_suffix}"
  temp_outfile2="${temp_outfile%.gz}"

  if [[ $temp_infile == *.gz ]]; then

    echo 'cp' "${temp_infile}" "${temp_outfile}"
    cp "${temp_infile}" "${temp_outfile}"
    echo 'bgzip -d' "${temp_outfile}"
    bgzip -d "${temp_outfile}"

    if [ "$cat_string" = "" ]; then
      cat_string="${temp_outfile2}"
    else
      cat_string="${cat_string} ${temp_outfile2}"
    fi

  else

    if [ "$cat_string" = "" ]; then
      cat_string="${temp_infile}"
    else
      cat_string="${cat_string} ${temp_infile}"
    fi

  fi

  echo ''
done

grep '^#' "${infile_for_header}" > "${infile_header}"

echo 'cat' $cat_string '| grep -v' '^#' '| cat' "${infile_header}" '- >' "${outfile}"
cat $cat_string | grep -v '^#' | cat "${infile_header}" - > "${outfile}"
echo ''

echo 'rm -rf' "${tmpdir}"'/'"${infile_prefix_basename}"'*'
rm -rf "${tmpdir}"/"${infile_prefix_basename}"*
echo ''

echo 'bgzip -f' "${outfile}"
bgzip -f "${outfile}"
echo ''

echo 'tabix -p vcf' "${outfile}"'.gz'
tabix -p vcf "${outfile}".gz
echo ''

echo ''
echo 'Finished! Output file is:' "${outfile}"'.gz'
echo ''

rm -rf $infile_header

touch "${done_file}"
rm -f "${lock_file}"
