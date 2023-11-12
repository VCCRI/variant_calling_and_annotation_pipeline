#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=2G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N samplesvcnv
#PBS -lstorage=gdata/abcd

set -euo pipefail

sample=$1
infile=$2
outdir=$3
outfile=$4
sw_and_refs=$5
large_homozygous_runs=$6
homozygous_runs=$7

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

outfile_basename=$(basename $outfile)
outfile_basename_prefix="${outfile_basename%.tsv}"
outfile_basename_prefix="${outfile_basename_prefix%.txt}"_dir

tmpdir="${PBS_JOBFS}"/scrap50/"${outfile_basename_prefix}"
mkdir -p "${tmpdir}"

# Where are the columns that contain the sample:
# ...INVBND2BDR5INS	INVBND2DEL	INVBND2OVERLAP	XSAMPLE_ONE	XSAMPLE_TWO	gene_CDSexons	gene_entireGene...
# ...INVBND2BDR5INS	INVBND2DEL	INVBND2OVERLAP	X4397	X4398	X974	VAF_in_INFO	VAF_4397	VAF_4398	VAF_974	gene_CDSexons	gene_entireGene	exon_containing_left_BND...
cols_string=$(head -n 1 $infile)
IFS=$'\t' read -r -a array <<< "$cols_string"
i=0
Xsample="X"$sample
vaf_sample="VAF_"$sample
sample_col=-1
vaf_col=-1
col_before_samples=-1
col_after_samples=-1
for this_col in "${array[@]}"; do
  i=$(( i + 1 ))
  IFS=$'.' read -r -a array2 <<< "$this_col"
  this_col_prefix="${array2[0]}"
  if [[ $this_col == $sample ]]; then
    if [[ "$sample_col" -eq -1 ]]; then
      sample_col=$i
    fi
  fi
  if [[ $this_col_prefix == $sample ]]; then
    if [[ "$sample_col" -eq -1 ]]; then
      sample_col=$i
    fi
  fi
  if [[ $this_col == $Xsample ]]; then
    sample_col=$i
  fi
  if [[ $this_col_prefix == $Xsample ]]; then
    sample_col=$i
  fi
  if [[ $this_col == $vaf_sample ]]; then
    vaf_col=$i
  fi
  if [[ $this_col_prefix == $vaf_sample ]]; then
    vaf_col=$i
  fi
  if [[ $this_col == "INVBND2OVERLAP" ]]; then
    col_before_samples=$i
  fi
  if [[ $this_col == "INV5" ]]; then
    col_before_samples=$i
  fi
  if [[ $this_col == "gene_CDSexons" ]]; then
    col_after_samples=$i
  fi
done

temp_sample_vars="${tmpdir}"/temp_sample_vars.txt
temp_outfile_large_hrun_flag="${tmpdir}"/temp_outfile_large_hrun_flag.txt
temp_outfile_hrun_flag="${tmpdir}"/temp_outfile_hrun_flag.txt

if [[ "$vaf_col" -eq -1 ]]; then
  vaf_col=sample_col
fi

echo 'cut -d$''\t' '-f"1-'$col_before_samples','$sample_col','$vaf_col','$col_after_samples'-"' $infile '>' $temp_sample_vars
cut -d$'\t' -f"1-$col_before_samples,$sample_col,$vaf_col,$col_after_samples-" $infile > $temp_sample_vars
echo ''

echo 'awk -f' $sw'/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file='"$large_homozygous_runs" '-v annot='"sample_BigHomozygRun" $temp_sample_vars '>' $temp_outfile_large_hrun_flag
awk -f $sw/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file="$large_homozygous_runs" -v annot="sample_BigHomozygRun" $temp_sample_vars > $temp_outfile_large_hrun_flag
echo ''

echo 'awk -f' $sw'/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file='"$homozygous_runs" '-v annot='"sample_HRun" $temp_outfile_large_hrun_flag '>' $outfile
awk -f $sw/victorchang_scripts/add_annot_region_file_annotation_to_tsv.awk -v annot_file="$homozygous_runs" -v annot="sample_HRun" $temp_outfile_large_hrun_flag > $outfile
echo ''

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

