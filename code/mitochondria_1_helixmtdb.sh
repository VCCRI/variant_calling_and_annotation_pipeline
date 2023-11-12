#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N helixmtdb
#PBS -lstorage=gdata/abcd

#infile=$1
#outdir=$2
#outfile=$3
#outfile_mito=$4
#sw_and_refs=$5

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

module load htslib/1.9

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

echo ''
echo 'Extract the mitochondrial variants'
echo ''

if [[ $infile == *.vcf.gz ]]; then
  echo 'tabix -h' $infile 'chrM >' $outfile_mito
  tabix -h $infile chrM > $outfile_mito
else
  echo 'grep -P' '^#|^chrM' $infile '>' $outfile_mito
  grep -P '^#|^chrM' $infile > $outfile_mito
fi

echo ''
echo "Start running annovar for" $infile
echo ''

echo $sw'/annovar/table_annovar_ed.pl' "$outfile_mito" ${humandb}'/ -vcfinput -buildver' $genome_version '\'
echo '  -out' "$outfile" '-remove \'
echo '  -protocol HelixMTdb_20200327 \'
echo '  -operation f -nastring . \'
echo '  -arg -time'

$sw/annovar/table_annovar_ed.pl "$outfile_mito" ${humandb}/ -vcfinput -buildver "$genome_version" \
  -out "$outfile" -remove \
  -protocol HelixMTdb_20200327 \
  -operation f -nastring . \
  -arg -time

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'outfile_mito' $outfile_mito
echo 'outfile' $outfile
echo ''
echo 'Finished!'
echo ''
