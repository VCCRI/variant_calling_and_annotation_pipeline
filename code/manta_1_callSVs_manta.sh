#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -l jobfs=1GB
#PBS -N manta
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

set -euo pipefail

#sample=$1
#infile=$2
#outdir=$3
#outdir2=$4
#outfile=$5
#sw_and_refs=$6

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

module load python3

. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outdir_manta="${outfile}"
outdir_manta_annotate="${outdir2}"
outfile_vcf_gz="${outdir_manta}"/results/variants/diploidSV.vcf.gz
outfile_copy_vcf_gz="${outdir_manta_annotate}"/"${sample}".manta_diploidSV.vcf.gz
outfile_copy_vcf="${outfile_copy_vcf_gz%.gz}"
outfile_with_inversions_prefix="${outdir_manta_annotate}"/"${sample}".manta_diploidSV_withInversions
outfile_with_inversions_vcf="${outfile_with_inversions_prefix}".vcf
outfile_with_inversions_tsv="${outfile_with_inversions_prefix}".tsv
outfile_with_inversions_tsv_rmvCols="${outfile_with_inversions_prefix}".rmvCols.tsv

echo ''
echo 'Run Manta to call SVs'
echo ''

$manta_configManta \
--bam "${infile}" \
--referenceFasta "${ref_fasta}" \
--runDir "${outdir_manta}"

"${outfile}"/runWorkflow.py

echo ''
echo 'Run Manta to call Inversions'
echo ''

cp "${outfile_vcf_gz}" "${outfile_copy_vcf_gz}"
gunzip -f "${outfile_copy_vcf_gz}"
$sw_convertInversion $sw_samtools "${ref_fasta}" "${outfile_copy_vcf}" > "${outfile_with_inversions_vcf}"

echo ''
echo 'Convert VCF to tab-delimited'
echo ''

python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_annovar_qual.py \
  -i "${outfile_with_inversions_vcf}" \
  -o "${outfile_with_inversions_tsv}"

echo ''
echo 'Rename columns to ensure all column names are unique'
echo ''

sed -e 's/\tEND\tSVLEN\tSVTYPE\t/\tSVEND\tSVLEN\tSVTYPE\t/' "${outfile_with_inversions_tsv}" > "${outfile_with_inversions_tsv_rmvCols}"

echo ''
echo 'Output file:' $outfile_with_inversions_tsv_rmvCols

echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

