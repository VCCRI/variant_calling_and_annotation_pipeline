#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l mem=8G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

# create a single text file with coverage output from picard wgs metrics
# from individual sample files

# input from /my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr/SAMPLE_ONE.markDup.setTags.bqsr.Picard_CollectWgsMetrics.txt

exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "Terminating.  $1" >&2
        exit 1
    fi
}

sw_and_refs=/my/directory/variant_calling_and_annotation_pipeline/code/where_are_softwares_and_references.sh
. "${sw_and_refs}"

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/bam_bqsr
outfile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".vcfQC_3_get_coverage.output.txt

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

:>$tmpfile
for infile in "${indir}"/*.markDup.setTags.bqsr.bam; do
  if [[ -L "$infile" ]]; then
    do_nothing=1 # This is a symlink to a file in a previous batch that will have already been QC-verified
  else
    sample=$(basename $infile)
    IFS='.' read -r -a array <<< "$sample"
    sample="${array[0]}"
    metrics_file="${indir}"/"${sample}".markDup.setTags.bqsr.Picard_CollectWgsMetrics.txt
    echo -e "${sample}\t${infile}\t${metrics_file}" >> "${tmpfile}"
  fi
done

manifest=$tmpfile
ls "$manifest" &>/dev/null
out_file=$outfile
exit_on_error "Unable to find manifest file: $manifest"
# store the subject id's in an array
a=( $( awk '{print $1}' $manifest ) )
metrics_files=( $( awk '{print $3}' $manifest ) )
# get the header from the first metrics file
s=${a[0]}
#f=$metrics_d/${s}$suf
f="${metrics_files[0]}"
ls "$f" &>/dev/null
exit_on_error "Unable to find metrics file: $f"
header=$( grep -E "^GENOME_TERRITORY" $f )
exit_on_error "Unable to get header for subject: $s"
idx=0
{
    printf "ID\t${header}\n"
    # print the header, and insert subject into first column
    for s in ${a[@]} ; do
        #f=$metrics_d/${s}$suf
        f="${metrics_files[$idx]}"
        coverage_lines="0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0"
        if [ -f $f ]; then
          ls "$f" &>/dev/null
          exit_on_error "Unable to find metrics file: $f"
          # the grep will return the header line + the data line for the coverage
          # Only want the data line
          coverage_lines=$( grep -A1 -E "^GENOME_TERRITORY" $f | tail -1 )
          exit_on_error "Unable to get coverage lines for subject: $s"
          # a the subject as the first column
          printf "$s\t${coverage_lines}\n"
        else
          printf "$s\t${coverage_lines}\n"
        fi
	(( idx=idx+1 ))
    done
} > $out_file

echo ''
echo 'outfile:' $outfile

echo ''
echo 'Finished!'
echo ''

#touch "${done_file}"
#rm -f "${lock_file}"



