#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=40G
#PBS -N spliceaiTensorflowResults
#PBS -lstorage=gdata/abcd

# Most sample chroms can process in walltime=5:00:00
# Some need more, eg. walltime=10:00:00

set -euo pipefail

sample=$1 # SAMPLE_ONE
infile=$2 # SAMPLE_ONE.annovar_clinvarDATE_spliceai.spliceogen.gnomad_filtered_splicing_variants.CHD_955_genes.needs_spliceai_tensorflow_scores.tensorflow_results.tsv
infile2=$3 # SAMPLE_ONE.annovar_clinvarDATE_spliceai.spliceogen.gnomad_filtered_splicing_variants.CHD_955_genes.needs_spliceai_tensorflow_scores.tsv
outdir=$4
outfile=$5
sw_and_refs=$6

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

module load bcftools
module load bedtools

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

##########

infile_basename=$(basename $infile)
outfile_basename=$(basename $outfile)
tmp_annot_file="${tmpdir}"/"${infile_basename}"
tmp_outfile="${tmpdir}"/"${outfile_basename}"

# Here is the infile header:
# CHROM	POS	END	ID	REF	ALT	QUAL	FILTER	INFO_SAMPLE_FOR_SPLICEAI	INFO_SpliceAI	INFO_SpliceAI_ALLELE	INFO_SpliceAI_SYMBOL	INFO_SpliceAI_DS_AG	INFO_SpliceAI_DS_ALINFO_SpliceAI_DS_DG	INFO_SpliceAI_DS_DL	INFO_SpliceAI_DP_AG	INFO_SpliceAI_DP_AL	INFO_SpliceAI_DP_DG	INFO_SpliceAI_DP_DL
# or
# CHROM	POS	END	ID	REF	ALT	QUAL	FILTER	SAMPLE_FOR_SPLICEAI	SpliceAI	SpliceAI_ALLELE	SpliceAI_SYMBOL	SpliceAI_DS_AG	SpliceAI_DS_ALSpliceAI_DS_DG	SpliceAI_DS_DL	SpliceAI_DP_AG	SpliceAI_DP_AL	SpliceAI_DP_DG	SpliceAI_DP_DL
# Convert it to:
# CHROM	POS	REF	ALT	SpliceAIresult	SpliceAIresult_ALLELE	SpliceAIresult_SYMBOL	SpliceAIresult_DS_AG	SpliceAIresult_DS_ALSpliceAIresult_DS_DG	SpliceAIresult_DS_DL	SpliceAIresult_DP_AG	SpliceAIresult_DP_AL	SpliceAIresult_DP_DG	SpliceAIresult_DP_DL

echo 'cut -d$''\t' '-f1-2,5-6,10-' $infile '| awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {gsub("INFO_SpliceAI","SpliceAIresult",$0); gsub("SpliceAI","SpliceAIresult",$0); gsub("SpliceAIresultresult","SpliceAIresult",$0); print $0} else {print $0}}' '>' $tmp_annot_file
cut -d$'\t' -f1-2,5-6,10- $infile | awk 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {gsub("INFO_SpliceAI","SpliceAIresult",$0); gsub("SpliceAI","SpliceAIresult",$0); gsub("SpliceAIresultresult","SpliceAIresult",$0); print $0} else {print $0}}' > $tmp_annot_file
echo ''

echo "Add the bed file of the spliceai tensorflow results to the tsv file as annotation columns"
echo ''

echo 'awk -f' $sw'/victorchang_scripts/add_annot_file_annotation_to_tsv.awk -v annot_file='"$tmp_annot_file" $infile2 '>' $tmp_outfile
awk -f $sw/victorchang_scripts/add_annot_file_annotation_to_tsv.awk -v annot_file="$tmp_annot_file" $infile2 > $tmp_outfile
echo ''

echo 'awk -v sample='"$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmp_outfile}" '>' "${outfile}"
awk -v sample="$sample" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {print $0, "Sample_Id"} else {print $0, sample}}' "${tmp_outfile}" > "${outfile}"
echo ''


touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo 'outfile:' $outfile
echo ''


