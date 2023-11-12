#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N pathogenic
#PBS -lstorage=gdata/abcd

sample=$1
infile=$2
outdir=$3
outfile=$4
sw_and_refs=$5

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

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# Use HelixMTdb for MAF values. Extract mitochondrial variants as rare when their allele frequency is ≤ 1%.
# From those, chose possibly_pathogenic variants as those that have at least one of the following:
# 1. MitImpact's APOGEE prediction for the variant is “P” for pathogenic; or
# 2. MitoTIP's prediction score for the variant excedes the recommended pathogenicity threshold of 12.66; or
# 3. The variant is marked as 'confirmed pathogenic' in MITOMAP.

apogee_col=-1
helix_hom_col=-1
helix_het_col=-1
mitotip_col=-1
mitomap_col=-1
extra1_col=-1
vaf_col=-1
gt_col=-1
cols_string=$(head -n 1 $infile)
IFS=$'\t' read -r -a array <<< "$cols_string"
i=0
for this_col in "${array[@]}"; do
  i=$(( i + 1 ))
  if [[ $this_col == "INFO_MitImpact_APOGEE" ]]; then
    apogee_col=$i
  fi
  if [[ $this_col == "MitImpact_APOGEE" ]]; then
    apogee_col=$i
  fi
  if [[ $this_col == "INFO_HelixMTdb_AF_hom" ]]; then
    helix_hom_col=$i
  fi
  if [[ $this_col == "HelixMTdb_AF_hom" ]]; then
    helix_hom_col=$i
  fi
  if [[ $this_col == "INFO_HelixMTdb_AF_het" ]]; then
    helix_het_col=$i
  fi
  if [[ $this_col == "HelixMTdb_AF_het" ]]; then
    helix_het_col=$i
  fi
  if [[ $this_col == "INFO_MitoTIP_Score" ]]; then
    mitotip_col=$i
  fi
  if [[ $this_col == "MitoTIP_Score" ]]; then
    mitotip_col=$i
  fi
  if [[ $this_col == "INFO_Mitomap_cal_aachange" ]]; then
    mitomap_col=$i
  fi
  if [[ $this_col == "Mitomap_cal_aachange" ]]; then
    mitomap_col=$i
  fi
  if [[ $this_col == "INFO_extra_mito_vars_consequence" ]]; then
    extra1_col=$i
  fi
  if [[ $this_col == "extra_mito_vars_consequence" ]]; then
    extra1_col=$i
  fi
  if [[ $this_col == "VAF" ]]; then
    vaf_col=$i
  fi
  if [[ $this_col == "FORMAT_GT" ]]; then
    gt_col=$i
  fi
done

echo 'awk -v sample='"$sample" '-v apogee_col='"$apogee_col" '-v helix_hom_col='"$helix_hom_col" '-v helix_het_col='"$helix_het_col" '-v mitotip_col='"$mitotip_col" '-v mitomap_col='"$mitomap_col" '-v extra1_col='"$extra1_col" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "sample"} else {if ( ($0 ~ /athogenic/) || ($apogee_col == "Pathogenic") || (($helix_hom_col <= 0.01) && ($helix_hom_col != ".") && ($vaf_col >= 0.65)) || (($helix_het_col <= 0.01) && ($helix_het_col != ".") && ($vaf_col <= 0.85)) || (($mitotip_col >= 12.66) && ($mitotip_col != ".")) || ($mitomap_col ~ /Pathogenic/) || ($mitomap_col ~ /confirmed_pathogenic/) || ($mitomap_col ~ /likely_pathogenic/) || (($extra1_col != ".")&&($extra1_col != ".")) ) {print $0}}}' $infile '>' $outfile

awk -v sample="$sample" -v apogee_col="$apogee_col" -v helix_hom_col="$helix_hom_col" -v helix_het_col="$helix_het_col" -v mitotip_col="$mitotip_col" -v mitomap_col="$mitomap_col" -v extra1_col="$extra1_col" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "sample"} else {if ( ($0 ~ /athogenic/) || ($apogee_col == "Pathogenic") || (($helix_hom_col <= 0.01) && ($helix_hom_col != ".") && ($vaf_col >= 0.65)) || (($helix_het_col <= 0.01) && ($helix_het_col != ".") && ($vaf_col <= 0.85)) || (($mitotip_col >= 12.66) && ($mitotip_col != ".")) || ($mitomap_col ~ /Pathogenic/) || ($mitomap_col ~ /confirmed_pathogenic/) || ($mitomap_col ~ /likely_pathogenic/) || (($extra1_col != ".")&&($extra1_col != ".")) ) {print $0, sample}}}' $infile > $outfile

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

