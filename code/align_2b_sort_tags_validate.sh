#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=8G
#PBS -l ncpus=2
#PBS -l jobfs=400G
#PBS -N sortTagValidate
#PBS -lstorage=gdata/abcd

JAVA_GC_THREADS=1
XMS=4
XMX=6

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
#sw_and_refs=where_are_softwares_and_references.sh
# sw_and_refs is passed in as a parameter
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# Mark duplicates already completed
# This program sorts, fixes tags, and validates BAM

tmpfile_sortSam="${tmpdir}"/"${sample}".gatk_sortSam_output.bam
infile_bytes=$( ls -l $infile | awk '{print $5}' )
# bytes in 200G - using base 10 system
BYTES_200G=$( python3 -c "print(200*10**9)" )
if [ $infile_bytes -gt $BYTES_200G ] ; then
    tmpfile_sortSam=${outdir}/${sample}_gatk_sortSam_output.bam
fi

echo ''
echo 'Run gatk SortSam on' $infile 'to produce' $tmpfile_sortSam
echo ''

#NOTE non-default values for
# COMPRESSION_LEVEL
# MAX_RECORDS_IN_RAM - just use default setting here, as high value will slow
# processing with smaller amounts of mem
# Also use of parallel garbage collection
${gatk_path} --java-options "-server -Xms${XMS}g -Xmx${XMX}g -Djava.io.tmpdir=$tmpdir -XX:+UseParallelGC -XX:ParallelGCThreads=$JAVA_GC_THREADS  " \
    SortSam \
    --INPUT "$infile" \
    --OUTPUT "$tmpfile_sortSam" \
    --SORT_ORDER "coordinate" \
    --TMP_DIR $tmpdir \
    --COMPRESSION_LEVEL 2 \
    --CREATE_INDEX false \
    --CREATE_MD5_FILE false

echo ''
echo 'Then run picard SetNmMdAndUqTags on' $tmpfile_sortSam 'to produce' $outfile
echo ''

# NOTE - leaving to default compression - as this is our final output file
${gatk_path} --java-options "-server -Xms${XMS}g -Xmx${XMX}g -Djava.io.tmpdir=$tmpdir -XX:+UseParallelGC -XX:ParallelGCThreads=$JAVA_GC_THREADS  " \
    SetNmMdAndUqTags \
    --INPUT "$tmpfile_sortSam" \
    --OUTPUT "$outfile" \
    --CREATE_INDEX true \
    --CREATE_MD5_FILE true \
    --REFERENCE_SEQUENCE ${ref_fasta}

# no errors - so remove temporary sort bam
# NOTE - this is particularly important if we are dealing with a large file that
# is using /g/data rather then jobfs
rm "$tmpfile_sortSam"

# Validate the BAM file

echo ''
echo 'Run picard ValidateSamFile on' $outfile 'to produce' "${outfile}".Picard_ValidateSamFile.txt
echo ''

# NOTE - in addition above parameter comments, using
# MAX_OPEN_TEMP_FILES=1000 as default is 8000 & it is thought our NCI ulimit is 1024
java -server -Xms${XMS}g -Xmx${XMX}g -Djava.io.tmpdir=$tmpdir -XX:+UseParallelGC -XX:ParallelGCThreads=$JAVA_GC_THREADS -jar $picard_jar \
    ValidateSamFile \
    I="$outfile" \
    TMP_DIR=$tmpdir \
    MAX_OPEN_TEMP_FILES=1000 \
    MODE=VERBOSE > "${outfile}".Picard_ValidateSamFile.txt

touch "${done_file}"
rm -f "${lock_file}"

echo 'Finished!'
echo ''
