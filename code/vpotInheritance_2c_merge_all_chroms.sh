#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=80G
#PBS -l ncpus=1
#PBS -l jobfs=2MB
#PBS -N catAllVpot
#PBS -lstorage=scratch/abcd+gdata/abcd
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

infile_prefix=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_vpot_ranking/"${batch}".gatk_vqsr_and_platypus_mnp.
infile_suffix=.annovar_clinvarDATE.vep_reformatted.extras.maf_0p01_and_scoring.vcf_final_output_file.adjust_scores.txt

outfile="${infile_prefix}AllChrom${infile_suffix}"

tmpdir="${currdir}"/tmp_for_vpotInheritance_2c_merge_all_chroms
mkdir -p "${tmpdir}"

mito="M"
mito22Y="22MY"
if [[ $genome_version == "hg19" ]]; then
  mito="MT"
  mito22Y="22MTY"
fi

cat_string=''
# 22/M/Y were merged to ensure that all samples have at least one variant in the file, otherwise vpot would have dropped samples and variants
for chrom in {1..21} "X" "$mito22Y" ; do
  new_file="${infile_prefix}chr${chrom}${infile_suffix}"
  if [[ $cat_string == '' ]]; then
    cat_string=$new_file
  else
    cat_string="${cat_string} ${new_file}"
  fi
done

tmphdr="${tmpdir}"/"${batch}".tmp_hdr.txt
vcf_anychrom="${infile_prefix}chr21${infile_suffix}"

echo 'grep -P' '^#|^Ranking' $vcf_anychrom '>' $tmphdr
grep -P '^#|^Ranking' $vcf_anychrom > $tmphdr
echo ''

echo 'cat' $cat_string '| grep -v' '^#' '| grep -v' '^Ranking' '| cat' $tmphdr '- >' $outfile
cat $cat_string | grep -v '^#' | grep -v '^Ranking' | cat $tmphdr - > $outfile
echo ''

echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

