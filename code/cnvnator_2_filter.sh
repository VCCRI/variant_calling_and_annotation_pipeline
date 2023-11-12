#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=02:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=2G
#PBS -N mantabedtools
#PBS -lstorage=gdata/abcd

set -euo pipefail

sample=$1
infile_tsv=$2
outdir=$3
outfile=$4
sw_and_refs=$5

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

module load R/3.6.1

random_number=$RANDOM
script_name=$(basename $0)
outfile_basename=$(basename $outfile)

# head SAMPLE1.markDup.setTags.bqsr.cnvnator.txt
# deletion        chr1:140001-179000      39000   0.708638      4.07422e-09     4.06424e-06     5.62111e-09     2.34676e-05     0.91181
# deletion        chr1:207001-258000      51000   0.0171175     3.12495e-12     0       	3.2525e-12      0       	0.673239
# duplication     chr1:260001-289000      29000   1.8726  	8.86881e-08     5.83485e-184    2.77497e-07     1.13536e-170    0.356906

echo ''
echo 'cnvnator_2_filter: filter cnvnator results so that all fields are < 0.05 (eval1_ttest, eval2_gauss_prob, eval3_ttest_midCNV, and eval4_gauss_prob_midCNV), and add heading'
echo ''

this_hdr="${tmpdir}"/temp."${sample}"."${script_name}".filtered."${random_number}".txt
tmp_outfile="${tmpdir}"/temp."${sample}"."${script_name}".tmp_outfile."${random_number}".txt

echo -e "#Chr\tstart\tend\tcnv_type\tcnv_size\tnormalized_read_depth\teval1_ttest\teval2_gauss_prob\teval3_ttest_midCNV\teval4_gauss_prob_midCNV\tfraction_reads_mapped_with_q0_quality" > "${this_hdr}"

Rscript $sw/victorchang_scripts/filter_cnvnator_output.R $infile_tsv $tmp_outfile
cat $tmp_outfile | awk 'BEGIN {FS="\t";OFS="\t"} {print $2,$1,$3,$4,$5,$6,$7,$8,$9}' | sed 's/-/:/' | awk 'BEGIN {FS=":";OFS="\t"} {print $1,$2,$3}' | cat $this_hdr - | cut -d$'\t' -f1-11 > $outfile

rm $this_hdr
rm $tmp_outfile

echo ''
echo 'cnvnator_2_filter output:' $outfile

echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"



