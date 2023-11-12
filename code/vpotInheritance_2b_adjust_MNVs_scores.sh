#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=10G
#PBS -l ncpus=1
#PBS -l jobfs=10GB
#PBS -N adjustVpotScoring
#PBS -lstorage=gdata/abcd

# vpot priority will be called after this job, to filter variants by MAF<0.01.
# When a sample in the joint-called vcf does not have any variants (eg. female samples have no Y chromosome variants),
# then vpot silently samples and drops other existing variants in other samples.
# To ensure this doesn't happen, let's merge the shorted chromosomes so that together all samples will have variants.

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
currdir=/my/directory/variant_calling_and_annotation_pipeline/code
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${currdir}"/tmp_for_vpotInheritance_2b_adjust_MNVs_scores
mkdir -p "${tmpdir}"

tmphdr="${tmpdir}"/"${batch}".tmp_hdr.txt
head -n 1 $infile > $tmphdr

# Sort by chrom+pos,
# assign the same highest ranking and score to overlapping variants such as MNVs and their singleton SNPs,
# then go back to being sorted by score and ranking.
echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $infile '| sort -k3,3V -k4,4V -k6,6V -k7,7V | \'
echo '  awk -f' $sw'/victorchang_scripts/adjust_vpot_score_for_MNVs.awk | sort -k2,2Vr -k1,1Vr -k3,3V -k4,4V -k6,6V -k7,7V | cat' $tmphdr '- >' $outfile
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $infile | sort -k3,3V -k4,4V -k6,6V -k7,7V | \
  awk -f $sw/victorchang_scripts/adjust_vpot_score_for_MNVs.awk | sort -k2,2Vr -k1,1Vr -k3,3V -k4,4V -k6,6V -k7,7V | cat $tmphdr - > $outfile
echo ''

echo 'outfile:' $outfile
echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"

