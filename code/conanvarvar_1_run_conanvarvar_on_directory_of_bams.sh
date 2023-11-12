#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l mem=128GB
#PBS -l ncpus=4
#PBS -N conanvarvar
#PBS -lstorage=gdata/abcd

set -euo pipefail

#in_bam_dir=$1
#outdir_sw_and_working_dir=$2
#outdir_results=$3
#sw_and_refs=$4

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

module load R/3.6.1
module unload intel-fc intel-cc
module load intel-compiler/2019.3.199
export R_LIBS_USER="${sw}"/conanvarvar/R_libraries_for_conanvarvar

cd "${outdir_sw_and_working_dir}"

# The default resolution bin size is 50kb
#         --binsize=50000 \
# High resolution bin size is 10kb
#         --binsize=10000 \

conanvarvar_reference=hg38
if [[ $genome_version == 'hg19' ]]; then
  conanvarvar_reference=hg19
fi

Rscript --vanilla "${outdir_sw_and_working_dir}"/ConanVarvar.R \
        --bamdir="${in_bam_dir}" \
        --reference="${conanvarvar_reference}" \
        --outdir="${outdir_results}" \
        --plotresults \
        --ncores=4 \
        --verbose

echo ''
echo 'Finished!'
echo ''

