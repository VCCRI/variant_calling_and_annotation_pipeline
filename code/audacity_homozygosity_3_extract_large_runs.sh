#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=20G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N extracts
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

set -euo pipefail

sample=$1
infile=$2
outdir=$3
outfile=$4
sw_and_refs=$5

#queue_file="${outfile}.queued"
#lock_file="${outfile}.lock"
#done_file="${outfile}.done"
#term_file="${outfile}.term"
#log_file="${outfile}.log"

#touch "${lock_file}"
#rm -f "${queue_file}"

tmpdir="${PBS_JOBFS}"
mkdir -p "${tmpdir}"

. "${sw_and_refs}"

module load R/4.0.0
module load htslib
module load bedtools
module load intel-compiler/2019.3.199


outfile_basename=$(basename $outfile)
outfile_png="${outfile%.txt}"
outfile_png="${outfile_png%.tsv}".png

# head /my/directory/variant_calling_and_annotation_pipeline/working_directory/audacity_homozygosity/SAMPLE_ONE_chr1/SAMPLE_ONE/SAMPLE_ONE_DIDOH3M2Regions.txt
# chr1	884537	925628	41091	0	84	5	NA
# chr1	929346	959193	29847	0	56	1	NA
# chr1	985268	996128	10860	0	27	0	NA

echo 'Rscript' $sw'/victorchang_scripts/identify_large_homozygous_runs_from_bed_file.R' $infile $outfile
Rscript $sw/victorchang_scripts/identify_large_homozygous_runs_from_bed_file.R $infile $outfile
echo ''

echo 'Rscript' $sw'/victorchang_scripts/plot_homozygous_runs_from_bed_file.R' $sample $genome_version $outfile $outfile_png
Rscript $sw/victorchang_scripts/plot_homozygous_runs_from_bed_file.R $sample $genome_version $outfile $outfile_png
echo ''

echo ''
echo 'outfile:' $outfile
echo 'outfile_png:' $outfile_png

echo ''

echo 'Finished!'
echo ''

#touch "${done_file}"
#rm -f "${lock_file}"



