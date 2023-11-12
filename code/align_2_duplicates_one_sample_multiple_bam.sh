#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l mem=100G
#PBS -l ncpus=4
#PBS -l jobfs=70G
#PBS -N markDup
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

. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# Mark duplicates

# for performance, create duplicate file on jobfs, then copy to /g/data
intermediate_outfile=${tmpdir}/"${sample}"dup.bam

echo ''
echo 'Run picard MarkDuplicates for' $infile_list 'to produce' $intermediate_outfile
echo ''

IFS=':' read -r -a array <<< "$infile_list"
input_string=""
for infile in "${array[@]}"; do
    if [[ "$input_string" == "" ]]; then
        input_string="INPUT="$infile
    else
        temp_string="${input_string} INPUT=${infile}"
        input_string=$temp_string
    fi
done

# NOTES - mictro
# - using 1 of 2 threads for Java parallel garbage collection ; see note above
# - using latest Picard - 2.22.0 as of run date - to include performance enhancements
# for MarkDuplicates
# - writing output to jobfs for speed - before copying to /g/data
# - MAX_RECORDS_IN_RAM=5000000 # this is 10x default value (only if use >=32G mem)
# - we don't use DT tag - do set CLEAR_DT to false
# - we don't need to add PG to every read record - so set ADD_PG_TAG_TO_READS=false
# - setting MAX_FILE_HANDLES=1000 (default 8000) ; because I believe our ulimit on NCI is 1024
# - compression level set to 2 (default 5) to save time (at expense of disk space)

java -server -Xms54g -Xmx60g -Djava.io.tmpdir=$tmpdir -XX:+UseParallelGC \
    -Dsamjdk.compression_level=2 -XX:ParallelGCThreads=1 -jar $picard_jar \
    MarkDuplicates \
    $input_string \
    OUTPUT="$outfile" \
    METRICS_FILE="${outfile}".MarkDuplicatesMetrics.txt \
    VALIDATION_STRINGENCY=SILENT \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    ASSUME_SORT_ORDER="queryname" \
    MAX_FILE_HANDLES=1000 \
    MAX_RECORDS_IN_RAM=5000000 \
    CLEAR_DT=false \
    ADD_PG_TAG_TO_READS=false \
    TMP_DIR=$tmpdir

# copy to /g/data
echo ''
echo "Copying output BAM file for MarkDuplicates to: $outfile"
cp "$intermediate_outfile" "$outfile"
echo ''

touch "${done_file}"
rm -f "${lock_file}"

echo 'Finished!'
echo ''
