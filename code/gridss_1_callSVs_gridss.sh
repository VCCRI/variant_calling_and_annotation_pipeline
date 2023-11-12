#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=60GB
#PBS -l jobfs=60GB
#PBS -l ncpus=16
#PBS -N gridss
#PBS -lstorage=gdata/abcd

set -euo pipefail

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

NCPUS=16

term_handler()
{
    rm -f "${lock_file}"
    touch "${term_file}"
    exit 1
}
trap 'term_handler' TERM

touch "${lock_file}"
rm -f "${queue_file}"

module load java/jdk-8.40
module load R/3.6.1
module load samtools/1.9
module load bwa/0.7.17

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

outfile_prefix="${outfile%.gz}"
outfile_prefix="${outfile_prefix%.vcf}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"
tmpdir_assembly="${outfile_prefix}".tmpdir_assembly
mkdir -p "${tmpdir_assembly}"
tmpdir_workingdir="${outfile_prefix}".tmpdir_workingdir
mkdir -p "${tmpdir_workingdir}"
logdir="${outfile_prefix}".logdir
mkdir -p "${logdir}"
log_file="${logdir}"/"${sample}".gridss.log

java \
    -ea -Xmx23G \
    -Dsamjdk.create_index=true \
    -Dsamjdk.use_async_io_read_samtools=true \
    -Dsamjdk.use_async_io_write_samtools=true \
    -Dsamjdk.use_async_io_write_tribble=true \
    -Dsamjdk.compression_level=1 \
    -cp $gridss_jar gridss.CallVariants \
    INPUT="${infile}" \
    INPUT_LABEL="${sample}" \
    OUTPUT="${outfile}" \
    ASSEMBLY="${tmpdir_assembly}"/assembly.bam \
    REFERENCE_SEQUENCE="${ref_fasta}" \
    THREADS=16 \
    TMP_DIR="${tmpdir}" \
    WORKING_DIR="${tmpdir_workingdir}" \
    BLACKLIST="${gridss_blacklist}" 2>&1 | tee -a "${log_file}"

echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

