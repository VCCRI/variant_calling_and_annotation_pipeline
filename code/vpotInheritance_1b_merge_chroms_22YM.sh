#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=50G
#PBS -l ncpus=1
#PBS -l jobfs=100G
#PBS -N merge22YM
#PBS -lstorage=scratch/abcd+gdata/bacd
#noPBS -m bea

# vpot priority will be called after this job, to filter variants by MAF<0.01.
# When a sample in the joint-called vcf does not have any variants (eg. female samples have no Y chromosome variants),
# then vpot silently samples and drops other existing variants in other samples.
# To ensure this doesn't happen, let's merge the shorted chromosomes so that together all samples will have variants.

set -euo pipefail

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

infile_prefix=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate/"${batch}".gatk_vqsr_and_platypus_mnp.chr
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.vcf
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate

outfile="${infile_prefix}22MY${infile_suffix}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

vcf_chrom22="${infile_prefix}22${infile_suffix}"
vcf_chromM="${infile_prefix}M${infile_suffix}"
vcf_chromY="${infile_prefix}Y${infile_suffix}"

if [[ $genome_version == "hg19" ]]; then
  vcf_chromM="${infile_prefix}MT${infile_suffix}"
fi

tmphdr="${tmpdir}"/"${batch}".tmp_hdr.txt

echo 'grep' '^#' $vcf_chrom22 '>' $tmphdr
grep '^#' $vcf_chrom22 > $tmphdr
echo ''

echo 'cat' $vcf_chrom22 $vcf_chromM $vcf_chromY '| grep -v' '^#' '| cat' $tmphdr '- >' $outfile
cat $vcf_chrom22 $vcf_chromM $vcf_chromY | grep -v '^#' | cat $tmphdr - > $outfile
echo ''

echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

