#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N converttotsv
#PBS -lstorage=gdata/abcd

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

module load python3/3.8.5

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outfile_prefix="${outfile%.tsv}"
outfile_prefix="${outfile_prefix%.txt}"
outfile_annovar="${outfile_prefix}"."${genome_version}"_multianno.vcf

echo $sw'/annovar/table_annovar_ed.pl' "$infile" ${humandb}'/ -vcfinput -buildver' $genome_version '\'
echo '  -out' "$outfile_prefix" '-remove \'
echo '  -protocol extramito \'
echo '  -operation f -nastring . \'
echo '  -arg -time'

$sw/annovar/table_annovar_ed.pl "$infile" ${humandb}/ -vcfinput -buildver $genome_version \
  -out "$outfile_prefix" -remove \
  -protocol extramito \
  -operation f -nastring . \
  -arg -time

echo ''
echo 'python3' $sw'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i' $outfile_annovar '-o' $outfile '-end_id YES'
python3 $sw/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i $outfile_annovar -o $outfile -end_id YES

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''
