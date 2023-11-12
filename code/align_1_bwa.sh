#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=32GB
#PBS -l ncpus=16
#PBS -l jobfs=1GB
#PBS -N bwa
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

NCPUS=16

module load samtools

. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

#sort_mem_per_thread=1792M # 28*1024/16
#java_opt="-server -Xms1g -Xmx8g" # for mem=32G case

#flowcell is from input
#lane is from input
ID=${flowcell}-${lane}
PU=${flowcell}.${lane}.${sample}
PL=ILLUMINA
SM=$sample
LB=$library

# $bwa and $ref_fasta are defined in $sw_and_refs (where_are_softwares_and_references.sh)
echo $bwa 'mem -t' $NCPUS '-R' "@RG\tID:$ID\tPU:$PU\tLB:$LB\tPL:$PL\tSM:$SM" "$ref_fasta" "$infile1" "$infile2"
echo 'followed by'
echo 'samtools view -1 - -o' $outfile
echo ''

$bwa mem \
    -t $NCPUS \
    -R "@RG\tID:$ID\tPU:$PU\tLB:$LB\tPL:$PL\tSM:$SM" \
    "$ref_fasta" "$infile1" "$infile2" | \
  samtools view -1 - -o $outfile

touch "${done_file}"
rm -f "${lock_file}"

echo 'Finished!'
echo ''
