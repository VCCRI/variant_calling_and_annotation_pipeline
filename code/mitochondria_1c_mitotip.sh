#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N mitotop
#PBS -lstorage=gdata/jb96+gdata/ra5+gdata/a32

# Most sample chroms can process in walltime=5:00:00
# Some need more, eg. walltime=10:00:00

infile=$1
outdir=$2
outfile=$3
sw_and_refs=$4

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

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo "Start running annovar for" $infile
echo ''

echo $sw'/annovar/table_annovar_ed.pl' "$infile" ${humandb}'/ -vcfinput -buildver' $genome_version '\'
echo '  -out' "$outfile" '-remove \'
echo '  -protocol mitotip_20201110 \'
echo '  -operation f -nastring . \'
echo '  -arg -time'

$sw/annovar/table_annovar_ed.pl "$infile" ${humandb}/ -vcfinput -buildver $genome_version \
  -out "$outfile" -remove \
  -protocol mitotip_20201110 \
  -operation f -nastring . \
  -arg -time

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'outfile:' "${outfile}"."${genome_version}"_multianno.vcf
echo ''
echo 'Finished!'
echo ''
