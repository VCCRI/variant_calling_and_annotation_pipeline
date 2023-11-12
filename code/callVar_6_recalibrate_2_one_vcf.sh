#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l mem=8G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -N ApplyVQSR
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

module load python3-as-python
module load java/jdk-8.40
module load R/3.6.1

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# snp/indel_filter_level values obtained from:
# https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.hg38.wgs.inputs.json
snp_filter_level=99.7
indel_filter_level=99.7

outfile_basename=$(basename $outfile)
tmp_outfile="${tmpdir}"/"${outfile_basename}"

echo ''
echo 'Run gatk indel ApplyVQSR' $infile 'using' $outfile_indel_recal 'to produce output file' $tmp_outfile
echo ''

${gatk_path} --java-options "-server -Xms1g -Xmx5g -Djava.io.tmpdir=$tmpdir" \
    ApplyVQSR \
    -O ${tmp_outfile} \
    -V ${infile} \
    --recal-file ${outfile_indel_recal} \
    --tranches-file ${outfile_indel_tranches} \
    --truth-sensitivity-filter-level ${indel_filter_level} \
    --create-output-variant-index true \
    -mode INDEL

echo ''
echo 'Run gatk indel ApplyVQSR' $tmp_outfile 'using' $outfile_snp_recal 'to produce output file' $outfile
echo ''

${gatk_path} --java-options "-server -Xms1g -Xmx5g -Djava.io.tmpdir=$tmpdir" \
    ApplyVQSR \
    -O ${outfile} \
    -V ${tmp_outfile} \
    --recal-file ${outfile_snp_recal} \
    --tranches-file ${outfile_snp_tranches} \
    --truth-sensitivity-filter-level ${snp_filter_level} \
    --create-output-variant-index true \
    -mode SNP

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

