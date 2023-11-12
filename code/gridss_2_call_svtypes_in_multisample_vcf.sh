#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=1:00:00
#PBS -l mem=2GB
#PBS -l jobfs=2GB
#PBS -l ncpus=1
#PBS -N svtypes
#PBS -lstorage=gdata/abcd

set -euo pipefail

# input variables
#sample=$1
#infile=$2
#outfile=$3
#sw_and_refs=$4

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

module load bedtools
module load python2-as-python

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

output_SV_vcf_file="${outfile}"
output_SV_tsv_file="${output_SV_vcf_file%.vcf}"
output_SV_tsv_file="${output_SV_tsv_file%.tsv}"
output_SV_tsv_file="${output_SV_tsv_file%.txt}"
output_SV_tsv_file="${output_SV_tsv_file}".tsv
output_SV_tsv_file_rmvCols="${output_SV_tsv_file}".rmvCols.tsv

GAP=3000
cohort_id=$sample

##################################################
# Gunzip the gridss output vcf file
##################################################

infile_gz="${tmpdir}"/"${sample}".gridss.vcf.gz
gridss_vcf_file="${tmpdir}"/"${sample}".gridss.vcf
echo ''
echo 'if necessary, gunzip' "${infile}" 'produce' "${gridss_vcf_file}"
echo ''
if [[ $infile =~ \.gz$ ]]; then
  cp "${infile}" "${infile_gz}"
  gunzip -f "${infile_gz}"
else
  cp "${infile}" "${gridss_vcf_file}"
fi

##################################################
# Filter out low quality SV calls and add TRANCHE2 quality value
##################################################

# For gridss 1.N.N output, add_tranche2_and_filter_gridss_multisample_vcf.py has no -vaf_version or has -vaf_version 2
# For gridss 2.N.N output, add_tranche2_and_filter_gridss_multisample_vcf.py has -vaf_version 2

echo ''
echo 'cat' "${gridss_vcf_file}" '| python' "${SVannotation_software_directory}"'/programs/add_tranche2_and_filter_gridss_multisample_vcf.py -max_BND_length_for_no_split_reads 100 -vaf_version 2 \'
echo '        >' "${tmpdir}/${sample}_gridss_BND_filtered1.vcf"
cat "${gridss_vcf_file}" | python "${SVannotation_software_directory}"/programs/add_tranche2_and_filter_gridss_multisample_vcf.py -max_BND_length_for_no_split_reads 100 -vaf_version 2 \
	> "${tmpdir}/${sample}_gridss_BND_filtered1.vcf"

echo ''
echo 'Make sure it is filtered by TRANCHE2'
echo ''
grep -P -v 'TRANCHE2=LOW' "${tmpdir}/${sample}_gridss_BND_filtered1.vcf" > "${tmpdir}/${sample}_gridss_BND_filtered.vcf" || true

##################################################
# Convert GRIDSS BND records into <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>
##################################################
echo 'For' ${sample} 'convert GRIDSS BND records into <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>'

echo 'sort input by GRIDSS EVENT'
grep '^#' "${tmpdir}/${sample}_gridss_BND_filtered.vcf" > "${tmpdir}/${sample}_gridss_hdr.vcf"
grep -v '^#' "${tmpdir}/${sample}_gridss_BND_filtered.vcf" | sort -t$'\t' -k 3,3 > "${tmpdir}/${sample}_gridss_nohdr_sortedByEvent.txt" || true
cat "${tmpdir}/${sample}_gridss_hdr.vcf" "${tmpdir}/${sample}_gridss_nohdr_sortedByEvent.txt" > "${tmpdir}/${sample}_gridss_sortedByEvent.vcf"

echo 'run program to identify <DEL>, <INS>, <INDEL>, <DUP:TANDEM>, <DUP:INS>'
echo 'cat' "${tmpdir}/${sample}_gridss_sortedByEvent.vcf" '\'
echo '        | python' "${SVannotation_software_directory}"'/programs/convert_multisample_VCF_BND_records_to_simple_structural_variant_VCF_records.py \'
echo '        >' "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf"
cat "${tmpdir}/${sample}_gridss_sortedByEvent.vcf" \
	| python "${SVannotation_software_directory}"/programs/convert_multisample_VCF_BND_records_to_simple_structural_variant_VCF_records.py \
	> "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf"

echo 'Remove SV variants having one or both breakends fall in a region of low-complexity or simple-repeats.'
echo 'They are probably mapped incorrectly and so the SV is false.'
echo 'These regions are the low-complexity or simple-repeat regions defined in UCSC Repeats file, not the entire Repeats file of mobile elements and other repeats.'
grep '^#' "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf" > "${tmpdir}/${cohort_id}_VCF_hdr.vcf"
grep -v '^#' "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_includes_LowComplexity.vcf" | sort -k1,1 -k2,2n > "${tmpdir}/${cohort_id}_VCF_body.vcf" || true

cat "${tmpdir}/${cohort_id}_VCF_body.vcf" | \
	python "${SVannotation_software_directory}"/programs/add_breakend_key_to_vcf_record.py -p LEFT > "${tmpdir}/${cohort_id}_add_key_LEFT.txt"
bedtools intersect -a "${tmpdir}/${cohort_id}_add_key_LEFT.txt" \
	-b "${regions_of_LowComplexity_SimpleRepeats}" -v | cut -d$'\t' -f4- > "${tmpdir}/${cohort_id}_add_key_LEFT_removed.vcf"
cat "${tmpdir}/${cohort_id}_add_key_LEFT_removed.vcf" | \
	python "${SVannotation_software_directory}"/programs/add_breakend_key_to_vcf_record.py -p RIGHT > "${tmpdir}/${cohort_id}_add_key_RIGHT.txt"
bedtools intersect -a "${tmpdir}/${cohort_id}_add_key_RIGHT.txt" \
	-b "${regions_of_LowComplexity_SimpleRepeats}" -v | cut -d$'\t' -f4- > "${tmpdir}/${cohort_id}_removed_low_complexity_body.vcf"
cat "${tmpdir}/${cohort_id}_VCF_hdr.vcf" "${tmpdir}/${cohort_id}_removed_low_complexity_body.vcf" > "${tmpdir}/${cohort_id}_removed_low_complexity.vcf"
rm "${tmpdir}/${cohort_id}_VCF_body.vcf"
rm "${tmpdir}/${cohort_id}_add_key_LEFT.txt" "${tmpdir}/${cohort_id}_add_key_LEFT_removed.vcf"
rm "${tmpdir}/${cohort_id}_add_key_RIGHT.txt" "${tmpdir}/${cohort_id}_removed_low_complexity_body.vcf"

cat "${tmpdir}/${cohort_id}_VCF_hdr.vcf" "${tmpdir}/${cohort_id}_removed_low_complexity.vcf" > "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP.vcf"

echo 'sort output by CHROM+POS'
grep -v '^#' "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP.vcf" > "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_nohdr.txt" || true
sort -t$'\t' -k1,1 -k2,2n "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP_nohdr.txt" > "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_nohdr_sorted.txt"

echo 'split output into the new INS_DEL_INDEL_DUP records and the remaining old BND records not used to call the INS_DEL_INDEL_DUP'
grep 'SVTYPE=BND' "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_nohdr_sorted.txt" > "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_nohdr_sorted.txt"
grep -v 'SVTYPE=BND' "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_nohdr_sorted.txt" > "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_nohdr_sorted.txt" || true
cat "${tmpdir}/${sample}_gridss_hdr.vcf" "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_nohdr_sorted.txt" > "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP.vcf"

echo 'identify new VCF headers added by INS_DEL_INDEL_DUP program'
grep '^#' "${tmpdir}/${sample}_gridss_sortedByEvent_BND_INS_DEL_INDEL_DUP.vcf" > "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_hdr.vcf"
sort "${tmpdir}/${sample}_gridss_hdr.vcf" | uniq > "${tmpdir}/${sample}_gridss_hdr_sorted.txt"
sort "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_hdr.vcf" | uniq > "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_hdr_sorted.txt"
comm -23 "${tmpdir}/${sample}_gridss_BND_INS_DEL_INDEL_DUP_hdr_sorted.txt" "${tmpdir}/${sample}_gridss_hdr_sorted.txt" > "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_new_hdrs_only.txt"

##################################################
# Convert GRIDSS BND records into <INV> and other inversion records
##################################################
echo 'For' ${sample} 'convert GRIDSS BND records into <INV> and other inversion records'

echo 'sort input by special key needed to put BND records of one inversion next to each other, in preparation for identifying inversions'
python "${SVannotation_software_directory}"/programs/add_key_to_VCF_to_find_inversions.py \
	-i "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP.vcf" \
	-o "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key.txt"
sort -k1,1 -k2,2n -k3,3 -k4,4n "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key.txt" > "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key_sorted.txt"

echo 'run program to identify inversions, but do not yet identify circular-inversions'
echo 'python' "${SVannotation_software_directory}"'/programs/identify_inversions_in_multisample_VCF_with_key.py \'
echo '        -i' "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key_sorted.txt" '\'
echo '        -o' "${tmpdir}/${sample}_gridss_inversions_only.vcf" '\'
echo '        -u' "${tmpdir}/${sample}_gridss_BND_not_used_for_inversions.vcf" '\'
echo '        -gap' "${GAP}"
python "${SVannotation_software_directory}"/programs/identify_inversions_in_multisample_VCF_with_key.py \
	-i "${tmpdir}/${sample}_gridss_BND_not_used_for_INS_DEL_INDEL_DUP_with_key_sorted.txt" \
	-o "${tmpdir}/${sample}_gridss_inversions_only.vcf" \
	-u "${tmpdir}/${sample}_gridss_BND_not_used_for_inversions.vcf" \
	-gap "${GAP}"

echo 'sort input by special key needed to put BND records of one inversion next to each other, in preparation for identifying circular-inversions'
python "${SVannotation_software_directory}"/programs/add_key_to_VCF_to_find_inversions.py \
	-i "${tmpdir}/${sample}_gridss_BND_not_used_for_inversions.vcf" \
	-o "${tmpdir}/${sample}_gridss_unused_breakend_records_with_key.txt"
sort -k1,1 -k2,2n -k3,3 -k4,4n "${tmpdir}/${sample}_gridss_unused_breakend_records_with_key.txt" > "${tmpdir}/${sample}_gridss_unused_breakend_records_with_key_sorted.txt"

echo 'run program to identify circular-inversions'
echo 'python' "${SVannotation_software_directory}"'/programs/identify_inversions_in_multisample_VCF_with_key.py \'
echo '        -i' "${tmpdir}/${sample}_gridss_unused_breakend_records_with_key_sorted.txt" '\'
echo '        -o' "${tmpdir}/${sample}_gridss_circular_inversions_too.vcf" '\'
echo '        -circular CIRCULAR \'
echo '        -u' "${tmpdir}/${sample}_gridss_BND_not_used_for_circular_inversions.vcf" '\'
echo '        -gap' "${GAP}"
python "${SVannotation_software_directory}"/programs/identify_inversions_in_multisample_VCF_with_key.py \
	-i "${tmpdir}/${sample}_gridss_unused_breakend_records_with_key_sorted.txt" \
	-o "${tmpdir}/${sample}_gridss_circular_inversions_too.vcf" \
	-circular CIRCULAR \
	-u "${tmpdir}/${sample}_gridss_BND_not_used_for_circular_inversions.vcf" \
	-gap "${GAP}"

echo 'identify new VCF headers added by inversion program'
echo 'the circular_inversions output has same VCF headers as inversions output plus a few more VCF headers'
grep '^#' "${tmpdir}/${sample}_gridss_circular_inversions_too.vcf" > "${tmpdir}/${sample}_gridss_circular_inversions_too_hdr.vcf"
sort "${tmpdir}/${sample}_gridss_circular_inversions_too_hdr.vcf" > "${tmpdir}/${sample}_gridss_circular_inversions_too_hdr_sorted.txt"
comm -23 "${tmpdir}/${sample}_gridss_circular_inversions_too_hdr_sorted.txt" "${tmpdir}/${sample}_gridss_hdr_sorted.txt" > "${tmpdir}/${sample}_gridss_inversions_new_hdrs_only.txt"

echo 'join INS_DEL_INDEL_DUP, inversions and circular-inversions output, and also include remaining unsed BND records, sort it all by CHROM+POS, add headers ready for input to genomic-region-annotation'
grep -v '^#' "${tmpdir}/${sample}_gridss_inversions_only.vcf" > "${tmpdir}/${sample}_gridss_inversions_only_nohdr.txt" || true
grep -v '^#' "${tmpdir}/${sample}_gridss_circular_inversions_too.vcf" > "${tmpdir}/${sample}_gridss_circular_inversions_too_nohdr.txt" || true
grep -v '^#' "${tmpdir}/${sample}_gridss_BND_not_used_for_circular_inversions.vcf" > "${tmpdir}/${sample}_gridss_BND_not_used_for_circular_inversions_nohdr.txt" || true
cat "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_nohdr_sorted.txt" "${tmpdir}/${sample}_gridss_inversions_only_nohdr.txt" "${tmpdir}/${sample}_gridss_circular_inversions_too_nohdr.txt" "${tmpdir}/${sample}_gridss_BND_not_used_for_circular_inversions_nohdr.txt" \
	| sort -k1,1 -k2,2n > "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_INV_BND_nohdr_sorted.txt"
grep '^#CHROM' "${tmpdir}/${sample}_gridss_hdr.vcf" > "${tmpdir}/${sample}_gridss_hdr_CHROM_only.txt"
grep -v '^#CHROM' "${tmpdir}/${sample}_gridss_hdr.vcf" > "${tmpdir}/${sample}_gridss_hdr_no_CHROM.txt" || true
cat "${tmpdir}/${sample}_gridss_hdr_no_CHROM.txt" "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_new_hdrs_only.txt" "${tmpdir}/${sample}_gridss_inversions_new_hdrs_only.txt" "${tmpdir}/${sample}_gridss_hdr_CHROM_only.txt" \
	"${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_INV_BND_nohdr_sorted.txt" > "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_INV_BND_got_GT_of_dot.vcf"

echo 'python' "${SVannotation_software_directory}"'/programs/convert_sample_GT_in_VCF_file.py -i' "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_INV_BND_got_GT_of_dot.vcf" '-o' "${output_SV_vcf_file}"' -a ALL'
python "${SVannotation_software_directory}"/programs/convert_sample_GT_in_VCF_file.py -i "${tmpdir}/${sample}_gridss_INS_DEL_INDEL_DUP_INV_BND_got_GT_of_dot.vcf" -o "${output_SV_vcf_file}" -a ALL

echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_annovar_qual.py -i' "${output_SV_vcf_file}" '-o' "${output_SV_tsv_file}"
python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_annovar_qual.py -i "${output_SV_vcf_file}" -o "${output_SV_tsv_file}"

##################################################
# Remove some of the GRIDSS output fields, and make sure remaining column names are unique
##################################################
echo 'For' ${sample} 'remove some of the GRIDSS output fields'

# #Chr	Start	End	Ref	Alt	Qual	AS	ASQ	ASRP	ASSR	BA	BAQ	BEID	BQ	BSC	BSCQ	BUM	BUMQ	CAS	CASQ	CIEND_1	CIEND_2	CIPOS_1	CIPOS_2	CIRPOS_1	CIRPOS_2	CQ	END	EVENT	HOMLEN	HOMSEQ	IC	IHOMPOS_1	IHOMPOS_2	IMPRECISE	IQ	PARID	RAS	RASQ	REF	REFPAIR	RP	RPQ	RSI	SC	SELF	SI	SR	SRQ	SVLEN	SVTYPE	TRANCHE	TRANCHE2	BNDVAF	CIEND_1	CIEND_2	CIPOS_1	CIPOS_2	INSSEQ	VAF	gridssAS	gridssASQ	gridssASRP	gridssASRR	gridssBA	gridssBAQ	gridssBEID	gridssBQ	gridssBSC	gridssBSCQ	gridssBUM	gridssBUMQ	gridssCAS	gridssCASQ	gridssCIEND	gridssCIPOS	gridssCIRPOSgridssCQ	gridssFILTER	gridssHOMLEN	gridssHOMSEQ	gridssIC	gridssIHOMPOS	gridssIQ	gridssQUAL	gridssRAS	gridssRASQ	gridssREF	gridssREFPAIR	gridssRP	gridssRPQ	gridssRSI	gridssSI	gridssSR	gridssSRQ	gridssVAF	INVBND1BDR3INS	INVBND1BDR5INS	INVBND1DEL	INVBND1OVERLAP	INVBND2BDR3INS	INVBND2BDR5INS	INVBND2DEL	INVBND2OVERLAP	sample_id
# cut -d$'\t' -f1-6,28,50- | sed -e 's/\tEND\tSVLEN\tSVTYPE\t/\tSVEND\tSVLEN\tSVTYPE\t/'
# #Chr	Start	End	Ref	Alt	Qual	SVEND	SVLEN	SVTYPE	TRANCHE	TRANCHE2	BNDVAF	CIEND_1	CIEND_2	CIPOS_1	CIPOS_2	INSSEQ	VAF	gridssAS	gridssASQ	gridssASRP	gridssASRR	gridssBA	gridssBAQ	gridssBEID	gridssBQ	gridssBSC	gridssBSCQ	gridssBUM	gridssBUMQ	gridssCAS	gridssCASQ	gridssCIEND	gridssCIPOS	gridssCIRPOS	gridssCQ	gridssFILTER	gridssHOMLEN	gridssHOMSEQ	gridssIC	gridssIHOMPOS	gridssIQ	gridssQUAL	gridssRAS	gridssRASQ	gridssREF	gridssREFPAIR	gridssRP	gridssRPQ	gridssRSI	gridssSI	gridssSR	gridssSRQ	gridssVAF	INVBND1BDR3INS	INVBND1BDR5INS	INVBND1DEL	INVBND1OVERLAP	INVBND2BDR3INS	INVBND2BDR5INS	INVBND2DEL	INVBND2OVERLAP	sample_id

echo 'cut -d$\t -f1-6,28,50-' "${output_SV_tsv_file}" '| sed -e s/\tEND\tSVLEN\tSVTYPE\t/\tSVEND\tSVLEN\tSVTYPE\t/ >' "${output_SV_tsv_file_rmvCols}"
cut -d$'\t' -f1-6,28,50- "${output_SV_tsv_file}" | sed -e 's/\tEND\tSVLEN\tSVTYPE\t/\tSVEND\tSVLEN\tSVTYPE\t/' > "${output_SV_tsv_file_rmvCols}"



echo 'Finished!'

touch "${done_file}"
rm -f "${lock_file}"

