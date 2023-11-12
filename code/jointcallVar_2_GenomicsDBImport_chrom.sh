#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=40:00:00
#PBS -l mem=30G
#PBS -l ncpus=4
#PBS -l jobfs=1G
#PBS -N GenomicsDBImport
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

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
tmpdir2="${PBS_JOBFS}"/tmp2
mkdir -p "${tmpdir}"
mkdir -p "${tmpdir2}"

echo ''
echo 'Run gatk GenomicsDBImport on chromosome' $chrom 'of' $infile_list 'to produce output genomicsdb-workspace-path directory' $outfile
echo ''

IFS=':' read -r -a array <<< "$infile_list"
input_string=""
for infile in "${array[@]}"; do
    if [[ "$input_string" == "" ]]; then
        input_string="-V "$infile
    else
        temp_string="${input_string} -V ${infile}"
        input_string=$temp_string
    fi
done

# This job runs GATK GenomicsDBImport to create a genomicsdb-workspace-path directory.
# The outdir is the genomicsdb-workspace-path directory
# and needs to be non-existent or empty directory before GATK GenomicsDBImport runs,
# so this job first deletes it if it exists already.
echo ''
echo 'Delete the directory' $outfile 'if it exists, before running GATK GenomicsDBImport to create it as a genomicsdb-workspace-path.'
echo ''
if [[ -f $outfile ]]; then
    echo 'The arguement passed as the destination for creating a genomicsdb-workspace-path directory with GATK GenomicsDBImport is in fact an existing file and thus I do not dare delete it and instead abort this job:' $outfile
    exit 1
fi
rm -rf "${outfile}"

echo ${gatk_path} '--java-options "-server -Xms1g -Xmx8g -Djava.io.tmpdir=$tmpdir" \'
echo '    GenomicsDBImport \'
echo '    '$input_string '\'
echo '    -L '${chrom} '\'
echo '    --genomicsdb-workspace-path' $outfile '\'
echo '    --tmp-dir $tmpdir2 \'
echo '    --reader-threads 4'

${gatk_path} --java-options "-server -Xms1g -Xmx8g -Djava.io.tmpdir=$tmpdir" \
    GenomicsDBImport \
    $input_string \
    -L ${chrom} \
    --genomicsdb-workspace-path $outfile \
    --tmp-dir $tmpdir2 \
    --reader-threads 4

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''
