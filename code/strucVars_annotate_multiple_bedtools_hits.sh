#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=64G
#PBS -l ncpus=1
#PBS -l jobfs=400G
#PBS -N strucVarsAnnotate
#PBS -lstorage=gdata/abcd

set -euo pipefail

#sample=$1
#infile_tsv=$2
#outdir=$3
#outfile=$4
#outprefix=$5
#sw_and_refs=$6

# This script uses bedtools to find the intersection between each row of the input data and the reference tables.
# Bedtools produces multiple output rows for a given input row when it overlaps with more than one entry in the reference table.
# Often the reference table value is the same value for all those overlaps. Eg. the reference table contains multiple transcript for the same region for the same gene.
# This script then collapses those multiple rows to one row for each different value it overlaps with.
# It does this by sorting, then running a script to fill in the unique value(s) of the last column(s) added by the bedtools intersect.
# Please note that this script uses linux sort on the first 5 columns instead of bedtools.
# Rows having the same beginning key are then collapsed to one row.
# For multiple rows of the following records where columns 1 to 3 are the same:
# chr17   77522483        81970118        T       <DUP:TANDEM>
# chr17   77522483        77522536        T       TATACACATATATATATATATATATATATACACACATATATATATATATATATAC
# bedtools sort does not sort on the fifth column.
# The consequences of this is that multiple rows remains, which can exponentially amplify the number of rows for each subsequent reference table comparison.

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

touch "${lock_file}"
rm -f "${queue_file}"

module load python
module load bedtools
module load R/3.6.1
module unload intel-fc intel-cc
module load intel-compiler/2019.3.199


# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outprefix_sorted="${outprefix}".sorted.tsv
outprefix_CDSexons="${outprefix}".CDSx.tsv
outprefix_CDSexons_collapsed="${outprefix}".CDSx.collapsed.tsv
outprefix_genes="${outprefix}".CDSx.Gene.tsv
outprefix_genes_collapsed="${outprefix}".CDSx.Gene.collapsed.tsv
outprefix_bndInsideExon="${outprefix}".CDSx.Gene.bndX.tsv
outprefix_gencode_CDSexons="${outprefix}".CDSx.Gene.bndX.gCDSx.tsv
outprefix_gencode_CDSexons_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.collapsed.tsv
outprefix_gencode_genes="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.tsv
outprefix_gencode_genes_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.collapsed.tsv
outprefix_gencode_bndInsideExon="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.tsv
outprefix_EHRFr104="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.tsv
outprefix_EHRFr104_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.collapsed.tsv
outprefix_gnomadsv="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.tsv
outprefix_gnomadsv_AF="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gnomadSV_AF.tsv
outprefix_gnomadsv_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.collapsed.tsv
outprefix_dgv="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.tsv
outprefix_dgv_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.collapsed.tsv
outprefix_segdup="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.tsv
outprefix_dbscSNV_largeSV="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.large_dbscSNV.tsv
outprefix_dbscSNV="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.tsv
outprefix_dbscSNV_1lastCol="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.1lastCol.tsv
outprefix_dbscSNV_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.collapsed.tsv
outprefix_FANTOM5_TSS="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.tsv
outprefix_FANTOM5_TSS_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.Dup.DGV.dbsc.TSS.collapsed.tsv
outprefix_pext_artery_aorta="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.px1.tsv
outprefix_pext_artery_aorta_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.Dup.DGV.dbsc.TSS.px1.collapsed.tsv
outprefix_pext_artery_coronary="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.px2.tsv
outprefix_pext_artery_coronary_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.Dup.DGV.dbsc.TSS.px2.collapsed.tsv
outprefix_pext_heart_atrial_appendage="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.px3.tsv
outprefix_pext_heart_atrial_appendage_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.Dup.DGV.dbsc.TSS.px3.collapsed.tsv
outprefix_pext_heart_left_ventricle="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.px4.tsv
outprefix_pext_heart_left_ventricle_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.Dup.DGV.dbsc.TSS.px4.collapsed.tsv
outprefix_pext_genes="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.tsv
outprefix_pext_genes_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.Dup.DGV.dbsc.TSS.pext.collapsed.tsv
outprefix_cleanIntronDels="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.tsv
outprefix_cleanIntronDels_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.collapsed.tsv
outprefix_intronTooShort="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.tsv
outprefix_intronTooShort_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.collapsed.tsv
outprefix_gencode_intronTooShort="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.tsv
outprefix_gencode_intronTooShort_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.collapsed.tsv
outprefix_hotspot_exons="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.hot.tsv
outprefix_hotspot_exons_collapsed="${outprefix}".CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.hot.collapsed.tsv

random_number=$RANDOM
script_name=$(basename $0)
outfile_basename=$(basename $outfile)



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: sort the input'
echo ''

this_input="${infile_tsv}"
echo 'sed s/^#Chr/0chrom/' "${this_input}" '| sed s/^#CHROM/0chrom/ | sed s/^CHROM/0chrom/ | sed s/^chrom/0chrom/ | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_sorted}"
sed 's/^#Chr/0chrom/' "${this_input}" | sed 's/^#CHROM/0chrom/' | sed 's/^CHROM/0chrom/' | sed 's/^chrom/0chrom/' | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_sorted}"
echo ''
this_input="${outprefix_sorted}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_refseq_cdsexons'
echo ''

# head $ucsc_refseq_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".CDSexons."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgene_CDSexons\tgene_CDSexons_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgene_CDSexons\tgene_CDSexons_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_refseq_cdsexons}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_CDSexons}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" | bedtools intersect -wao -a - -b "${ucsc_refseq_cdsexons}" | \
  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat "${this_hdr}" - | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_CDSexons}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It might be needed for any duplicate exons, otherwise same gene might appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_CDSexons}" '>' "${outprefix_CDSexons_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_CDSexons}" > "${outprefix_CDSexons_collapsed}"
echo ''
this_input="${outprefix_CDSexons_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_refseq_genes'
echo ''

# head $ucsc_refseq_genes
# chr1	11873	14409	DDX11L1	+
# chr1	14361	29370	WASH7P	-
# chr1	17368	17436	MIR6859-1	-
# chr1	29925	31295	MIR1302-2HG	+
# chr1	30365	30503	MIR1302-2	+

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".genes."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 23
cut1_plus_4=$((cut1+4)) # 27
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgene_entireGene\tgene_entireGene_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgene_entireGene\tgene_entireGene_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_refseq_genes}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_genes}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" | bedtools intersect -wao -a - -b "${ucsc_refseq_genes}" | \
  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat "${this_hdr}" - | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_genes}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_genes}" '>' "${outprefix_genes_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_genes}" > "${outprefix_genes_collapsed}"
echo ''
this_input="${outprefix_genes_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do left BND fall inside an exon?'
echo ''

# head $ucsc_refseq_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

bndInsideExon_left_bnd="${tmpdir}"/temp."${sample}".bndInsideExon_left_bnd."${outfile_basename}".bed
bndInsideExon_left_bnd_intersect="${tmpdir}"/temp."${sample}".bndInsideExon_left_bnd_intersect."${outfile_basename}".bed
bndInsideExon_left_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".bndInsideExon_left_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tleft_BND_inside_exon\texon_containing_left_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tleft_BND_inside_exon\texon_containing_left_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${bndInsideExon_left_bnd}"
awk 'BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" | cat "${this_hdr}" - > "${bndInsideExon_left_bnd}" || true
echo 'bedtools intersect -wao -a' "${bndInsideExon_left_bnd}" '-b' "${ucsc_refseq_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- >' "${bndInsideExon_left_bnd_intersect}"
bedtools intersect -wao -a "${bndInsideExon_left_bnd}" -b "${ucsc_refseq_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - > "${bndInsideExon_left_bnd_intersect}"
echo ''
rm $bndInsideExon_left_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${bndInsideExon_left_bnd_intersect}" '>' "${bndInsideExon_left_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${bndInsideExon_left_bnd_intersect}" > "${bndInsideExon_left_bnd_intersect_collapsed}"
echo ''
rm $bndInsideExon_left_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do right BND fall inside an exon?'
echo ''

bndInsideExon_right_bnd="${tmpdir}"/temp."${sample}".bndInsideExon_right_bnd."${outfile_basename}".bed
bndInsideExon_right_bnd_intersect="${tmpdir}"/temp."${sample}".bndInsideExon_right_bnd_intersect."${outfile_basename}".bed
bndInsideExon_right_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".bndInsideExon_right_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tright_BND_inside_exon\texon_containing_right_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tright_BND_inside_exon\texon_containing_right_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${bndInsideExon_right_bnd}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" | cat "${this_hdr}" - > "${bndInsideExon_right_bnd}" || true
echo 'bedtools intersect -wao -a' "${bndInsideExon_right_bnd}" '-b' "${ucsc_refseq_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${bndInsideExon_right_bnd_intersect}"
bedtools intersect -wao -a "${bndInsideExon_right_bnd}" -b "${ucsc_refseq_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${bndInsideExon_right_bnd_intersect}"
echo ''
rm $bndInsideExon_right_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${bndInsideExon_right_bnd_intersect}" '>' "${bndInsideExon_right_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${bndInsideExon_right_bnd_intersect}" > "${bndInsideExon_right_bnd_intersect_collapsed}"
echo ''
rm $bndInsideExon_right_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: merge BND_inside_exon left BND and right BND into one file'
echo ''

echo 'awk -v left_bnd_file='"${bndInsideExon_left_bnd_intersect_collapsed}" '\'
echo '    -v right_bnd_file='"${bndInsideExon_right_bnd_intersect_collapsed}"' -f' $sw'/victorchang_scripts/combine_bnd_files_key5cols.awk \'
echo '    '"${this_input}" '>' "${outprefix_bndInsideExon}"
awk -v left_bnd_file="${bndInsideExon_left_bnd_intersect_collapsed}" \
    -v right_bnd_file="${bndInsideExon_right_bnd_intersect_collapsed}" -f $sw/victorchang_scripts/combine_bnd_files_key5cols.awk \
    "${this_input}" > "${outprefix_bndInsideExon}"
echo ''
rm $bndInsideExon_left_bnd_intersect_collapsed
rm $bndInsideExon_right_bnd_intersect_collapsed
this_input="${outprefix_bndInsideExon}"



if [[ "$ucsc_gencode_cdsexons" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_gencode_cdsexons'
echo ''

# head $ucsc_gencode_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gencode_CDSexons."${random_number}".txt
tmpfile11="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile11."${random_number}".txt
tmpfile12="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile12."${random_number}".txt
tmpfile13="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile13."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgencode_gene_CDSexons\tgencode_gene_CDSexons_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgencode_gene_CDSexons\tgencode_gene_CDSexons_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_gencode_cdsexons}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_gencode_CDSexons}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" > $tmpfile11
bedtools intersect -wao -a $tmpfile11 -b "${ucsc_gencode_cdsexons}" > $tmpfile12
rm $tmpfile11
sort $tmpfile12 | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 > $tmpfile13
rm $tmpfile12
cat "${this_hdr}" $tmpfile13 | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_gencode_CDSexons}"
rm "${this_hdr}"
rm $tmpfile13
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It might be needed for any duplicate exons, otherwise same gene might appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_gencode_CDSexons}" '>' "${outprefix_gencode_CDSexons_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_gencode_CDSexons}" > "${outprefix_gencode_CDSexons_collapsed}"
echo ''
this_input="${outprefix_gencode_CDSexons_collapsed}"

fi



if [[ "$ucsc_gencode_genes" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with ucsc_gencode_genes'
echo ''

# head $ucsc_refseq_genes
# chr1	11873	14409	DDX11L1	+
# chr1	14361	29370	WASH7P	-
# chr1	17368	17436	MIR6859-1	-
# chr1	29925	31295	MIR1302-2HG	+
# chr1	30365	30503	MIR1302-2	+

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".genes."${random_number}".txt
tmpfile21="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile21."${random_number}".txt
tmpfile22="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile22."${random_number}".txt
tmpfile23="${tmpdir}"/temp."${sample}"."${script_name}".tmpfile23."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 23
cut1_plus_4=$((cut1+4)) # 27
echo 'head -n 1' "${this_input}" '| sed -e s/$/\t.\t.\t.\tgencode_gene_entireGene\tgencode_gene_entireGene_strand/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed -e 's/$/\t.\t.\t.\tgencode_gene_entireGene\tgencode_gene_entireGene_strand/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${ucsc_gencode_genes}" '| \'
echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- | cut -d$\t -f"1-'$cut1','$cut1_plus_4 '>' "${outprefix_gencode_genes}"
echo 'awk' $this_input 'to' $tmpfile21
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" > $tmpfile21
bedtools intersect -wao -a $tmpfile21 -b "${ucsc_gencode_genes}" > $tmpfile22
rm $tmpfile21
sort $tmpfile22 | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | uniq > $tmpfile23
rm $tmpfile22
cat "${this_hdr}" $tmpfile23 | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" > "${outprefix_gencode_genes}"
rm "${this_hdr}"
rm $tmpfile23
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_gencode_genes}" '>' "${outprefix_gencode_genes_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_gencode_genes}" > "${outprefix_gencode_genes_collapsed}"
echo ''
this_input="${outprefix_gencode_genes_collapsed}"

fi



if [[ "$ucsc_gencode_cdsexons" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do left BND fall inside a gencode exon?'
echo ''

# head $ucsc_refseq_cdsexons
# chr1	67096251	67096321	C1orf141	-
# chr1	67103237	67103382	C1orf141	-
# chr1	67111576	67111644	C1orf141	-

gencode_bndInsideExon_left_bnd="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_left_bnd."${outfile_basename}".bed
gencode_bndInsideExon_left_bnd_intersect="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_left_bnd_intersect."${outfile_basename}".bed
gencode_bndInsideExon_left_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_left_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gencode_bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tgencode_left_BND_inside_exon\tgencode_exon_containing_left_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tgencode_left_BND_inside_exon\tgencode_exon_containing_left_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${gencode_bndInsideExon_left_bnd}"
awk 'BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" | cat "${this_hdr}" - > "${gencode_bndInsideExon_left_bnd}" || true
echo 'bedtools intersect -wao -a' "${gencode_bndInsideExon_left_bnd}" '-b' "${ucsc_gencode_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- >' "${gencode_bndInsideExon_left_bnd_intersect}"
bedtools intersect -wao -a "${gencode_bndInsideExon_left_bnd}" -b "${ucsc_gencode_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - > "${gencode_bndInsideExon_left_bnd_intersect}"
echo ''
rm $gencode_bndInsideExon_left_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${gencode_bndInsideExon_left_bnd_intersect}" '>' "${gencode_bndInsideExon_left_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${gencode_bndInsideExon_left_bnd_intersect}" > "${gencode_bndInsideExon_left_bnd_intersect_collapsed}"
echo ''
rm $gencode_bndInsideExon_left_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do right BND fall inside a gencode exon?'
echo ''

gencode_bndInsideExon_right_bnd="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_right_bnd."${outfile_basename}".bed
gencode_bndInsideExon_right_bnd_intersect="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_right_bnd_intersect."${outfile_basename}".bed
gencode_bndInsideExon_right_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".gencode_bndInsideExon_right_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gencode_bndInsideExon_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tgencode_right_BND_inside_exon\tgencode_exon_containing_right_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tgencode_right_BND_inside_exon\tgencode_exon_containing_right_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${gencode_bndInsideExon_right_bnd}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" | cat "${this_hdr}" - > "${gencode_bndInsideExon_right_bnd}" || true
echo 'bedtools intersect -wao -a' "${gencode_bndInsideExon_right_bnd}" '-b' "${ucsc_gencode_cdsexons}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${gencode_bndInsideExon_right_bnd_intersect}"
bedtools intersect -wao -a "${gencode_bndInsideExon_right_bnd}" -b "${ucsc_gencode_cdsexons}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $8}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${gencode_bndInsideExon_right_bnd_intersect}"
echo ''
rm $gencode_bndInsideExon_right_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes longer to run. 
# It is needed for the duplicate exons, otherwise same gene will appear multiple times in output column.
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values.
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${gencode_bndInsideExon_right_bnd_intersect}" '>' "${gencode_bndInsideExon_right_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${gencode_bndInsideExon_right_bnd_intersect}" > "${gencode_bndInsideExon_right_bnd_intersect_collapsed}"
echo ''
rm $gencode_bndInsideExon_right_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: merge BND_inside_exon left BND and right BND into one file'
echo ''

echo 'awk -v left_bnd_file='"${gencode_bndInsideExon_left_bnd_intersect_collapsed}" '\'
echo '    -v right_bnd_file='"${gencode_bndInsideExon_right_bnd_intersect_collapsed}"' -f' $sw'/victorchang_scripts/combine_bnd_files_key5cols.awk \'
echo '    '"${this_input}" '>' "${outprefix_gencode_bndInsideExon}"
awk -v left_bnd_file="${gencode_bndInsideExon_left_bnd_intersect_collapsed}" \
    -v right_bnd_file="${gencode_bndInsideExon_right_bnd_intersect_collapsed}" -f $sw/victorchang_scripts/combine_bnd_files_key5cols.awk \
    "${this_input}" > "${outprefix_gencode_bndInsideExon}"
echo ''
rm $gencode_bndInsideExon_left_bnd_intersect_collapsed
rm $gencode_bndInsideExon_right_bnd_intersect_collapsed
this_input="${outprefix_gencode_bndInsideExon}"

fi




echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with EHRFr104 (Ensemble Human Regulatory Features)'
echo ''

# head $refdata_EHRFr104
# chrom	start	end	feature_type
# chr2	113379801	113380200	Promoter Flanking Region
# chr18	32661402	32662400	Promoter Flanking Region
# chr3	41288801	41289200	CTCF Binding Site

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".EHRFr104_header."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 17
cut1_plus_4=$((cut1+4)) # 21
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tensemble_human_regulatory_features/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tensemble_human_regulatory_features/' > "${this_hdr}"

echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" '| bedtools intersect -wao -a - -b' "${refdata_EHRFr104}" '| cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_EHRFr104}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' "${this_input}" | bedtools intersect -wao -a - -b "${refdata_EHRFr104}" | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_EHRFr104}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_EHRFr104}" '>' "${outprefix_EHRFr104_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_EHRFr104}" > "${outprefix_EHRFr104_collapsed}"
echo ''
this_input="${outprefix_EHRFr104_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with gnomAD-SV data, must have 50% for both sample-SV and gnomad-SV'
echo ''

# head $refdata_gnomadsv
# chr1	10000	20000	gnomAD_v2_DUP_1_1	10000	DUP	20175	0.9395080208778381	21474	1217	9479	41	0.955695986747742	NA	NA	NA	BAF,RD	depth	NA	NA	NA	NA	NA	NA	NA	NA	OR4F5	NA	NA	True
# chr1	10642	10643	gnomAD_v2_BND_1_1	-1	BND	17	0.0008430000161752105	20178	7	5	10077	0.00116999994497746	NA	NA	NA	PE,SR	manta	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	False
# chr1	14500	43500	gnomAD_v2_DUP_1_2	29000	DUP	259	0.25	1036	161	49	308	0.417452991008759	NA	NA	NA	BAF,RDdepth	NA	NA	NA	NA	NA	NA	NA	NA	OR4F5	NA	NA	True

# 70% overlap is chosen because the Database of Genomic Variants http://dgv.tcag.ca/dgv/app/about?ref=GRCh37/hg19 uses 70% for merging SVs/CNVs as being the same variant.

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".gnomadsv_header."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 18
cut1_plus_1=$((cut1+1)) # 19
cut1_plus_2=$((cut1+2)) # 20
cut1_plus_3=$((cut1+3)) # 21
cut1_plus_4=$((cut1+4)) # 22
cut1_plus_5=$((cut1+5)) # 23
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\tgnomadsv_CHROM\tgnomadsv_START\tgnomadsv_END\tgnomadsv_NAME\tgnomadsv_SVLEN\tgnomadsv_SVTYPE\tgnomadsv_AC\tgnomadsv_AF\tgnomadsv_AN\tgnomadsv_N_HET\tgnomadsv_N_HOMALT\tgnomadsv_N_HOMREF\tgnomadsv_POPMAX_AF\tgnomadsv_CPX_INTERVALS\tgnomadsv_CPX_TYPE\tgnomadsv_STRANDS\tgnomadsv_EVIDENCE\tgnomadsv_ALGORITHMS\tgnomadsv_SOURCE\tgnomadsv_PROTEIN_CODING__COPY_GAIN\tgnomadsv_PROTEIN_CODING__DUP_LOF\tgnomadsv_PROTEIN_CODING__DUP_PARTIAL\tgnomadsv_PROTEIN_CODING__INTRONIC\tgnomadsv_PROTEIN_CODING__INV_SPAN\tgnomadsv_PROTEIN_CODING__LOF\tgnomadsv_PROTEIN_CODING__MSV_EXON_OVR\tgnomadsv_PROTEIN_CODING__NEAREST_TSS\tgnomadsv_PROTEIN_CODING__PROMOTER\tgnomadsv_PROTEIN_CODING__UTR\tgnomadsv_PROTEIN_CODING__INTERGENIC\tgnomadsv_column// >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\tgnomadsv_CHROM\tgnomadsv_START\tgnomadsv_END\tgnomadsv_NAME\tgnomadsv_SVLEN\tgnomadsv_SVTYPE\tgnomadsv_AC\tgnomadsv_AF\tgnomadsv_AN\tgnomadsv_N_HET\tgnomadsv_N_HOMALT\tgnomadsv_N_HOMREF\tgnomadsv_POPMAX_AF\tgnomadsv_CPX_INTERVALS\tgnomadsv_CPX_TYPE\tgnomadsv_STRANDS\tgnomadsv_EVIDENCE\tgnomadsv_ALGORITHMS\tgnomadsv_SOURCE\tgnomadsv_PROTEIN_CODING__COPY_GAIN\tgnomadsv_PROTEIN_CODING__DUP_LOF\tgnomadsv_PROTEIN_CODING__DUP_PARTIAL\tgnomadsv_PROTEIN_CODING__INTRONIC\tgnomadsv_PROTEIN_CODING__INV_SPAN\tgnomadsv_PROTEIN_CODING__LOF\tgnomadsv_PROTEIN_CODING__MSV_EXON_OVR\tgnomadsv_PROTEIN_CODING__NEAREST_TSS\tgnomadsv_PROTEIN_CODING__PROMOTER\tgnomadsv_PROTEIN_CODING__UTR\tgnomadsv_PROTEIN_CODING__INTERGENIC\tgnomadsv_column/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${refdata_gnomadsv}" '-f 0.50 -r | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-" | sed s/^0chrom/chrom/ >' "${outprefix_gnomadsv}"
bedtools intersect -wao -a "${this_input}" -b "${refdata_gnomadsv}" -f 0.50 -r | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-" | sed 's/^0chrom/chrom/' > "${outprefix_gnomadsv}"
rm "${this_hdr}"
echo 'cut -d$\t -f"1-'$cut1','$cut1_plus_3','$cut1_plus_5'"' "${outprefix_gnomadsv}" '| \'
echo '  awk -v cut1_plus_1='"$cut1_plus_1"' -v cut1_plus_2='"$cut1_plus_2"' BEGIN {FS="\t";OFS=""} {print $0 "\t", '$cut1_plus_1' "_" '$cut1_plus_2'} | \'
echo '  sed s/\._\./\./g | cut -d$\t -f"1-'$cut1','$cut1_plus_3'" | sed s/gnomadsv_SVTYPE_gnomadsv_AF/gnomadsv_SVTYPE_AF/ \'
echo '  sed s/^0chrom/0chrom/ | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_gnomadsv_AF}"
cut -d$'\t' -f"1-$cut1,$cut1_plus_3,$cut1_plus_5" "${outprefix_gnomadsv}" | \
  awk -v cut1_plus_1="$cut1_plus_1" -v cut1_plus_2="$cut1_plus_2" 'BEGIN {FS="\t";OFS=""} {print $0 "\t", $cut1_plus_1 "_" $cut1_plus_2}' | \
  sed 's/\._\./\./g' | cut -d$'\t' -f"1-$cut1,$cut1_plus_3" | sed 's/gnomadsv_SVTYPE_gnomadsv_AF/gnomadsv_SVTYPE_AF/' | \
  sed 's/^chrom/0chrom/' | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_gnomadsv_AF}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_gnomadsv_AF}" '>' "${outprefix_gnomadsv_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_gnomadsv_AF}" > "${outprefix_gnomadsv_collapsed}"
echo ''
this_input="${outprefix_gnomadsv_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with DGV data, must have 50% for both sample-SV and gnomad-SV'
echo ''

# head $refdata_DGVgoldStdVar_r20160516
# chr1	61723	297968	CNV_Loss_Coe2014,Vogler2010
# chr1	70007	88084	CNV_Gain_Perry2008,Sudmant2013
# chr1	86778	91967	CNV_Gain_Perry2008,Conrad2009,Sudmant2013

# 70% overlap is chosen because the Database of Genomic Variants http://dgv.tcag.ca/dgv/app/about?ref=GRCh37/hg19 uses 70% for merging SVs/CNVs as being the same variant.

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".DGV_header."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 19
cut1_plus_4=$((cut1+4)) # 23
echo 'head -n 1 "${this_input}" | sed s/^chrom/0chrom/ | sed -e s/$/\tdgv_chrom\tdgv_start\tdgv_end\tdgv_cnvtype_study/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\tdgv_chrom\tdgv_start\tdgv_end\tdgv_cnvtype_study/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${refdata_DGVgoldStdVar_r20160516}" '-f 0.50 -r | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_dgv}"
bedtools intersect -wao -a "${this_input}" -b "${refdata_DGVgoldStdVar_r20160516}" -f 0.50 -r | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_dgv}"
rm "${this_hdr}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_dgv}" '>' "${outprefix_dgv_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_dgv}" > "${outprefix_dgv_collapsed}"
echo ''
this_input="${outprefix_dgv_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do left BND fall in a segmental duplication region?'
echo ''

# head $refdata_segmental_duplications
# chr1	10169	37148	chr1:180723	0	+	chr1	180723	207666	26943	1	1000	N/A	N/A	N/A	N/A	align_both/0014/both0071547	27025	30	128	26897	26628	269	164	105	0.989998884633974	0.988895903739741	0.0100683956884972	0.0100742688854825
# chr1	180723	207666	chr1:10169	0	+	chr1	10169	37148	26979	1	1000	N/A	N/A	N/A	N/A	align_both/0014/both0071547	27025	30	128	26897	26628	269	164	105	0.989998884633974	0.988895903739741	0.0100683956884972	0.0100742688854825
# chr1	88000	121417	chr1:265774	0	+	chr1	265774	297956	32182	2	1000	N/A	N/A	N/A	N/A	align_both/0014/both0071548	33449	25	1299	32150	31941	209	133	76	0.993499222395023	0.992727272727273	0.00652911487615718	0.00653207256400941

left_bnd="${tmpdir}"/temp."${sample}".left_bnd."${outfile_basename}".bed
left_bnd_intersect="${tmpdir}"/temp."${sample}".left_bnd_intersect."${outfile_basename}".bed
left_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".left_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".segdup_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tleft_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_left_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tleft_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_left_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS=""} {if (NR>1) print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${left_bnd}"
awk 'BEGIN {FS="\t";OFS=""} {if (NR>1) {print $1 "\t" $2 "\t" $2 "\t" $1 ":" $2 "-" $3 ":" $4 ":" $5}}' "${this_input}" | cat "${this_hdr}" - > "${left_bnd}" || true
echo 'bedtools intersect -wao -a' "${left_bnd}" '-b' "${refdata_segmental_duplications}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${left_bnd_intersect}"
bedtools intersect -wao -a "${left_bnd}" -b "${refdata_segmental_duplications}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${left_bnd_intersect}"
echo ''
rm $left_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${left_bnd_intersect}" '>' "${left_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${left_bnd_intersect}" > "${left_bnd_intersect_collapsed}"
echo ''
rm $left_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: do right BND fall in a segmental duplication region?'
echo ''

right_bnd="${tmpdir}"/temp."${sample}".right_bnd."${outfile_basename}".bed
right_bnd_intersect="${tmpdir}"/temp."${sample}".right_bnd_intersect."${outfile_basename}".bed
right_bnd_intersect_collapsed="${tmpdir}"/temp."${sample}".right_bnd_intersect_collapsed."${outfile_basename}".bed
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".segdup_bedtools_header."${random_number}".txt
echo 'echo -e "chrom\tstart\tend\tright_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_right_BND" >' "${this_hdr}"
echo -e "chrom\tstart\tend\tright_BND_in_segmental_duplication_region\tsegmental_duplication_region_hit_by_right_BND" > "${this_hdr}"
echo 'awk BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" '| cat' "${this_hdr}" '- >' "${right_bnd}"
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $1 "\t" $3 "\t" $3 "\t" $1 ":" $2 "-" $3 ":" $4 ":"$5}}' "${this_input}" | cat "${this_hdr}" - > "${right_bnd}" || true
echo 'bedtools intersect -wao -a' "${right_bnd}" '-b' "${refdata_segmental_duplications}" '| \'
echo '  awk BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7} | cat' "${this_hdr}" '- | sed s/:-1--1//g >' "${right_bnd_intersect}"
bedtools intersect -wao -a "${right_bnd}" -b "${refdata_segmental_duplications}" | \
  awk 'BEGIN {FS="\t";OFS=""} {print $1, "\t", $2, "\t", $3, "\t", $4, "\t", $5 ":" $6 "-" $7}' | cat "${this_hdr}" - | sed 's/:-1--1//g' > "${right_bnd_intersect}"
echo ''
rm $right_bnd

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${right_bnd_intersect}" '>' "${right_bnd_intersect_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${right_bnd_intersect}" > "${right_bnd_intersect_collapsed}"
echo ''
rm $right_bnd_intersect

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: merge segmental duplication left BND and right BND into one file'
echo ''

echo 'awk -v left_bnd_file='"${left_bnd_intersect_collapsed}" '\'
echo '    -v right_bnd_file='"${right_bnd_intersect_collapsed}"' -f' $sw'/victorchang_scripts/combine_bnd_files_key5cols.awk \'
echo '    '"${this_input}" '>' "${outprefix_segdup}"
awk -v left_bnd_file="${left_bnd_intersect_collapsed}" \
    -v right_bnd_file="${right_bnd_intersect_collapsed}" -f $sw/victorchang_scripts/combine_bnd_files_key5cols.awk \
    "${this_input}" > "${outprefix_segdup}"
echo ''
rm $left_bnd_intersect_collapsed
rm $right_bnd_intersect_collapsed
this_input="${outprefix_segdup}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with dbscSNV splice sites'
echo ''

# head $refdata_dbscSNV
# #Chr	Start	End	Ref	Alt	dbscSNV_ada	dbscSNV_rf
# chr1	925799	925799	A	G	0.987833212454089	0.616
# chr1	925799	925799	A	T	0.976226938154799	0.578
# chr1	925800	925800	G	A	0.999648482686128	0.966

# Large SVs will hit a lot of dbscSNV entries and make an intermediate file that is too large to sort.
# We are not interested in the millions of splice site variants in large SVs anyway.
# So, for SVs greater than 1,000,000 bp, don't do the bedtools intersect with dbscSNV database.

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".dbscSNV_header."${random_number}".txt
empty_refdata_dbscSNV="${tmpdir}"/empty_refdata_dbscSNV."${sample}"."${script_name}".dbscSNV_header."${random_number}".txt
head -n 2 $refdata_dbscSNV > $empty_refdata_dbscSNV

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l) # 22
cut1_plus_1=$((cut1+1)) # 23
cut1_plus_2=$((cut1+2)) # 24
cut1_plus_3=$((cut1+3)) # 25
cut1_plus_4=$((cut1+4)) # 26
cut1_plus_5=$((cut1+5)) # 27
cut1_plus_6=$((cut1+6)) # 28
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\t.\t.\tdbscSNV_ada\tdbscSNV_rf\tbedtools_count/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\t.\t.\tdbscSNV_ada\tdbscSNV_rf\tbedtools_count/' > "${this_hdr}"

#bedtools intersect -wao -a "${this_input}" -b "${refdata_dbscSNV}" | cat "${this_hdr}" - | sed 's/^chrom/0chrom/' | \
#  cut -d$'\t' -f"1-$cut1,$cut1_plus_1-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_dbscSNV}"
#
echo 'awk BEGIN {FS="\t";OFS="\t"} {if ((NR!=1) && (($3-$2)>1000000)) {print $0}}' $this_input '| bedtools intersect -wao -a - -b' "${empty_refdata_dbscSNV}" '>' $outprefix_dbscSNV_largeSV
awk 'BEGIN {FS="\t";OFS="\t"} {if ((NR!=1) && (($3-$2)>1000000)) {print $0}}' $this_input | bedtools intersect -wao -a - -b "${empty_refdata_dbscSNV}" > $outprefix_dbscSNV_largeSV
#
echo 'awk BEGIN {FS="\t";OFS="\t"} {if ((NR==1) || (($3-$2)<=1000000)) {print $0}}' $this_input '| bedtools intersect -wao -a - -b' "${refdata_dbscSNV}" '| \'
echo '  cat' "${this_hdr}" '-' $outprefix_dbscSNV_largeSV '| sed s/^chrom/0chrom/ | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_1-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_dbscSNV}"
awk 'BEGIN {FS="\t";OFS="\t"} {if ((NR==1) || (($3-$2)<=1000000)) {print $0}}' $this_input | bedtools intersect -wao -a - -b "${refdata_dbscSNV}" | \
  cat "${this_hdr}" - $outprefix_dbscSNV_largeSV | sed 's/^chrom/0chrom/' | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_1-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_dbscSNV}"
#
rm "${this_hdr}" $empty_refdata_dbscSNV
echo ''
#
awk -v cut1_plus_1="$cut1_plus_1" -v cut1_plus_2="$cut1_plus_2" -v cut1_plus_4="$cut1_plus_4" -v cut1_plus_5="$cut1_plus_5" 'BEGIN {FS="\t";OFS=""} {print $0 "\t", $cut1_plus_1 ":" $cut1_plus_2 ":" $cut1_plus_4 ":" $cut1_plus_5}' "${outprefix_dbscSNV}" | sed 's/\.:-1:\.:./\./g' | cut -d$'\t' -f"1-$cut1,$cut1_plus_6" | sed 's/\t\.:\.:\.:\./\tdbscSNV/' > "${outprefix_dbscSNV_1lastCol}"
echo ''

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk' "${outprefix_dbscSNV_1lastCol}" '>' "${outprefix_dbscSNV_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique.awk "${outprefix_dbscSNV_1lastCol}" > "${outprefix_dbscSNV_collapsed}"
echo ''
this_input="${outprefix_dbscSNV_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with FANTOM5 CAGEr human heart transcription start sites'
echo ''

# head $refdata_FANTOM5_TSS
# #chrom	start	end	FANTOM5_heart_TSS_start	FANTOM5_heart_TSS_end	FANTOM5_heart_TSS_strand	FANTOM5_heart_TSS_tpm	FANTOM5_heart_TSS_gene
# chr1	630049	630150	630049	630150	+	0	MTND2P28
# chr1	632164	632199	632164	632199	+	17.62759588	MTCO1P12
# chr1	633533	633561	633533	633561	+	64.51133112	MTATP8P1

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".FANTOM5_TSS_header."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
cut1_plus_8=$((cut1+8))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tFANTOM5_heart_TSS_start/tFANTOM5_heart_TSS_end\tFANTOM5_heart_TSS_strand\tFANTOM5_heart_TSS_tpm\tFANTOM5_heart_TSS_gene/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tFANTOM5_heart_TSS_start\tFANTOM5_heart_TSS_end\tFANTOM5_heart_TSS_strand\tFANTOM5_heart_TSS_tpm\tFANTOM5_heart_TSS_gene/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${refdata_FANTOM5_TSS_5cols}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-$cut1_plus_8" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_FANTOM5_TSS}"
bedtools intersect -wao -a "${this_input}" -b "${refdata_FANTOM5_TSS_5cols}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-$cut1_plus_8" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_FANTOM5_TSS}" || true
rm "${this_hdr}"
echo ''

num_fantom5_cols=$(head -n 1 "${refdata_FANTOM5_TSS_5cols}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_FANTOM5_TSS}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_fantom5_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -v key=1,2,3,4,5 -v collapse='$list_columns_to_collapse '-f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk' "${outprefix_FANTOM5_TSS}" '>' "${outprefix_FANTOM5_TSS_collapsed}"
awk -v key=1,2,3,4,5 -v collapse=$list_columns_to_collapse -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk "${outprefix_FANTOM5_TSS}" > "${outprefix_FANTOM5_TSS_collapsed}"
echo ''
this_input="${outprefix_FANTOM5_TSS_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with PEXT transcript aware scores for artery_aorta'
echo ''

# #Chr	Start	End	PEXT_Artery_Aorta_exon_symbol	PEXT_Artery_Aorta_exon_avg_score
# chr1	135800	135810	AL627309.1	1
# chr1	137614	137623	AL627309.1	1
# chr1	138531	139309	AL627309.1	1
# chr1	924946	924956	SAMD11	0.256211
# chr1	925188	925197	SAMD11	0.0419255

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".pext_artery_aorta."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
cut1_plus_5=$((cut1+5))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tPEXT_Artery_Aorta_exon_symbol\tPEXT_Artery_Aorta_exon_avg_score/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tPEXT_Artery_Aorta_exon_symbol\tPEXT_Artery_Aorta_exon_avg_score/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${pext_artery_aorta}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_pext_artery_aorta}"
bedtools intersect -wao -a "${this_input}" -b "${pext_artery_aorta}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_pext_artery_aorta}" || true
rm "${this_hdr}"
echo ''

num_pext_cols=$(head -n 1 "${pext_artery_aorta}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_pext_artery_aorta}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_pext_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -v key=1,2,3,4,5 -v collapse='$list_columns_to_collapse '-f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk' "${outprefix_pext_artery_aorta}" '>' "${outprefix_pext_artery_aorta_collapsed}"
awk -v key=1,2,3,4,5 -v collapse=$list_columns_to_collapse -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk "${outprefix_pext_artery_aorta}" > "${outprefix_pext_artery_aorta_collapsed}"
echo ''
this_input="${outprefix_pext_artery_aorta_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with PEXT transcript aware scores for artery_coronary'
echo ''

# #Chr	Start	End	PEXT_Artery_Coronary_exon_symbol	PEXT_Artery_Coronary_exon_avg_score
# chr1	135800	135810	AL627309.1	1
# chr1	137614	137623	AL627309.1	1
# chr1	138531	139309	AL627309.1	1
# chr1	924946	924956	SAMD11	0.18956
# chr1	925798	925808	SAMD11	0.0714286

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".pext_artery_coronary."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
cut1_plus_5=$((cut1+5))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tPEXT_Artery_Coronary_exon_symbol\tPEXT_Artery_Coronary_exon_avg_score/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tPEXT_Artery_Coronary_exon_symbol\tPEXT_Artery_Coronary_exon_avg_score/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${pext_artery_coronary}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_pext_artery_coronary}"
bedtools intersect -wao -a "${this_input}" -b "${pext_artery_coronary}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_pext_artery_coronary}" || true
rm "${this_hdr}"
echo ''

num_pext_cols=$(head -n 1 "${pext_artery_coronary}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_pext_artery_coronary}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_pext_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -v key=1,2,3,4,5 -v collapse='$list_columns_to_collapse '-f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk' "${outprefix_pext_artery_coronary}" '>' "${outprefix_pext_artery_coronary_collapsed}"
awk -v key=1,2,3,4,5 -v collapse=$list_columns_to_collapse -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk "${outprefix_pext_artery_coronary}" > "${outprefix_pext_artery_coronary_collapsed}"
echo ''
this_input="${outprefix_pext_artery_coronary_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with PEXT transcript aware scores for heart_atrial_appendage'
echo ''

# #Chr	Start	End	PEXT_Heart_AtrialAppendage_exon_symbol	PEXT_Heart_AtrialAppendage_exon_avg_score
# chr1	135800	135810	AL627309.1	1
# chr1	137614	137623	AL627309.1	1
# chr1	138531	139309	AL627309.1	1
# chr1	924946	924956	SAMD11	0.229913
# chr1	925188	925197	SAMD11	0.00988875

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".pext_heart_atrial_appendage."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
cut1_plus_5=$((cut1+5))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tPEXT_Heart_AtrialAppendage_exon_symbol\tPEXT_Heart_AtrialAppendage_exon_avg_score/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tPEXT_Heart_AtrialAppendage_exon_symbol\tPEXT_Heart_AtrialAppendage_exon_avg_score/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${pext_heart_atrial_appendage}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_pext_heart_atrial_appendage}"
bedtools intersect -wao -a "${this_input}" -b "${pext_heart_atrial_appendage}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_pext_heart_atrial_appendage}" || true
rm "${this_hdr}"
echo ''

num_pext_cols=$(head -n 1 "${pext_heart_atrial_appendage}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_pext_heart_atrial_appendage}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_pext_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -v key=1,2,3,4,5 -v collapse='$list_columns_to_collapse '-f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk' "${outprefix_pext_heart_atrial_appendage}" '>' "${outprefix_pext_heart_atrial_appendage_collapsed}"
awk -v key=1,2,3,4,5 -v collapse=$list_columns_to_collapse -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk "${outprefix_pext_heart_atrial_appendage}" > "${outprefix_pext_heart_atrial_appendage_collapsed}"
echo ''
this_input="${outprefix_pext_heart_atrial_appendage_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with PEXT transcript aware scores for heart_left_ventricle'
echo ''

# #Chr	Start	End	PEXT_Heart_LeftVentricle_exon_symbol	PEXT_Heart_LeftVentricle_exon_avg_score
# chr1	135800	135810	AL627309.1	1
# chr1	137614	137623	AL627309.1	1
# chr1	138531	139309	AL627309.1	1
# chr1	930312	930344	SAMD11	0.888889
# chr1	931032	931097	SAMD11	0.888889

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".pext_heart_left_ventricle."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
cut1_plus_5=$((cut1+5))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tPEXT_Heart_LeftVentricle_exon_symbol\tPEXT_Heart_LeftVentricle_exon_avg_score/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tPEXT_Heart_LeftVentricle_exon_symbol\tPEXT_Heart_LeftVentricle_exon_avg_score/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${pext_heart_left_ventricle}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_pext_heart_left_ventricle}"
bedtools intersect -wao -a "${this_input}" -b "${pext_heart_left_ventricle}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-$cut1_plus_5" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_pext_heart_left_ventricle}" || true
rm "${this_hdr}"
echo ''

num_pext_cols=$(head -n 1 "${pext_heart_left_ventricle}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_pext_heart_left_ventricle}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_pext_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -v key=1,2,3,4,5 -v collapse='$list_columns_to_collapse '-f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk' "${outprefix_pext_heart_left_ventricle}" '>' "${outprefix_pext_heart_left_ventricle_collapsed}"
awk -v key=1,2,3,4,5 -v collapse=$list_columns_to_collapse -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk "${outprefix_pext_heart_left_ventricle}" > "${outprefix_pext_heart_left_ventricle_collapsed}"
echo ''
this_input="${outprefix_pext_heart_left_ventricle_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with PEXT transcript aware scores for genes (spanning the entire gene even though PEXT scores are for exons, so as to give an idea of how frequently this gene is expressed in this tissue)'
echo ''

# #Chr	Start	End	PEXT_gene_symbol	PEXT_gene_avg_Artery_Coronary	PEXT_gene_avg_Artery_Aorta	PEXT_gene_avg_Heart_AtrialAppendage	PEXT_gene_avg_Heart_LeftVentricle	PEXT_gene_avg_mean_proportion
# chr1	69091	70008	OR4F5	0	0	0	0	0
# chr1	135800	139309	AL627309.1	1	1	1	1	1
# chr1	450738	451676	OR4F29	0	0	0	0	1
# chr1	685716	686654	OR4F16	0	0	0	0	1
# chr1	803152	803757	AL669831.1	0	0	0	0	0
# chr1	882663	884603	AL645608.2	0	0	0	0	0
# chr1	924946	944153	SAMD11	0.641463	0.602555	0.626179	0.799824	0.694072
# chr1	925884	931065	AL645608.1	0	0	0	0	0
# chr1	944694	959240	NOC2L	0.762379	0.781371	0.759638	0.759249	0.731452

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".pext_genes."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
cut1_plus_9=$((cut1+9))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\tPEXT_gene_symbol\tPEXT_gene_avg_Artery_Coronary\tPEXT_gene_avg_Artery_Aorta\tPEXT_gene_avg_Heart_AtrialAppendage\tPEXT_gene_avg_Heart_LeftVentricle\tPEXT_gene_avg_mean_proportion/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\tPEXT_gene_symbol\tPEXT_gene_avg_Artery_Coronary\tPEXT_gene_avg_Artery_Aorta\tPEXT_gene_avg_Heart_AtrialAppendage\tPEXT_gene_avg_Heart_LeftVentricle\tPEXT_gene_avg_mean_proportion/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${pext_genes}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4-$cut1_plus_9" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_pext_genes}"
bedtools intersect -wao -a "${this_input}" -b "${pext_genes}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4-$cut1_plus_9" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_pext_genes}" || true
rm "${this_hdr}"
echo ''

num_pext_cols=$(head -n 1 "${pext_genes}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_pext_genes}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_pext_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -v key=1,2,3,4,5 -v collapse='$list_columns_to_collapse '-f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk' "${outprefix_pext_genes}" '>' "${outprefix_pext_genes_collapsed}"
awk -v key=1,2,3,4,5 -v collapse=$list_columns_to_collapse -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_by_column.awk "${outprefix_pext_genes}" > "${outprefix_pext_genes_collapsed}"
echo ''
this_input="${outprefix_pext_genes_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with introns and identify clean_intron_deletions'
echo ''

# head 
# chr1    67093605        67096250        C1orf141        -
# chr1    67096322        67103236        C1orf141        -
# chr1    33540699        33541128        CSMD2   -
# chr1    33541310        33542718        CSMD2   -

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".cleanIntronDels_header."${random_number}".txt
hdr_bit1=`head -n 1 "${this_input}" | sed 's/^chrom/0chrom/'`
hdr_bit2=$(echo -e "chrom\tstart\tend\tis_clean_intron_deletion")
echo -e "${hdr_bit1}\t${hdr_bit2}\toverlap" > "${this_hdr}"

chrom1=1
pos1=2
end1=3
alt1=5
chrom1=$(( $chrom1 + 3 ))
pos1=$(( $pos1 + 3 ))
end1=$(( $end1 + 3 ))
alt1=$(( $alt1 + 3 ))
num_cols1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
num_cols2=$(head -n 1 "${ucsc_refseq_introns}" | sed -e 's/\t/\n/g' | wc -l)
overlap=$(( $num_cols1 + 3 + $num_cols2 + 1 ))
chrom2=$(( overlap - 5 ))
pos2=$(( overlap - 4 ))
end2=$(( overlap - 3 ))
gene=$(( overlap - 2 ))

this_input_basename=$(basename $this_input)
tmpfile_cleanIntronDels="${tmpdir}"/"${this_input_basename}"
echo 'Create a temporary file of only the clean_intron_deletions for this sample.'
echo 'grep -v' 'SVTYPE=BND' "${this_input}" '| awk' 'BEGIN {FS="\t";OFS="\t"} {print $1, $2, $end1, $0}' '| awk' 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' '| bedtools intersect -a - -b' $ucsc_refseq_introns '-wao | \'
echo '  awk -v overlap='"$overlap" '-v pos1='"$pos1" '-v end1='"$end1" '-v alt1='"$alt1" '-v gene='"$gene" '-v chrom2='"$chrom2" '-v pos2='"$pos2" '-v end2='"$end2" 'function abs(v) {return v < 0 ? -v : v} BEGIN {FS="\t";OFS="\t"} {if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {print $1, $pos1, $end1, "is_clean_intron"}}' '>' $tmpfile_cleanIntronDels
grep -v 'SVTYPE=BND' "${this_input}" | awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $2, $3, $0}' | awk 'BEGIN {FS="\t";OFS="\t"} {if ($2 != $3) {print $0}}' | bedtools intersect -a - -b $ucsc_refseq_introns -wao | \
  awk -v overlap="$overlap" -v pos1="$pos1" -v end1="$end1" -v alt1="$alt1" -v gene="$gene" -v chrom2="$chrom2" -v pos2="$pos2" -v end2="$end2" 'function abs(v) {return v < 0 ? -v : v} BEGIN {FS="\t";OFS="\t"} {if ( ($overlap>0) && (abs($overlap-($end2-$pos2))<=5) && (abs(($end1-$pos1)-($end2-$pos2))<=5) && (abs($pos1-$pos2)<=5) && (abs($end1-$end2)<=5) ) {print $1, $pos1, $end1, "is_clean_intron"}}' > $tmpfile_cleanIntronDels || true
echo ''

num_clean_intron_deletions=$(wc -l "${tmpfile_cleanIntronDels}" | cut -d" " -f1)

if [[ $num_clean_intron_deletions -eq 0 ]]; then
  echo 'There are no clean_intron_deletions for this sample.'
  echo -e "${hdr_bit1}\tis_clean_intron_deletion" > "${this_hdr}"
  echo 'grep -v' '^chrom' "${this_input}" '| awk' 'BEGIN {FS="\t";OFS="\t"} {print $0, "."}' '| cat' "${this_hdr}" '- | sed s/^0chrom/chrom/ >' "${outprefix_cleanIntronDels}"
  grep -v '^chrom' "${this_input}" | awk 'BEGIN {FS="\t";OFS="\t"} {print $0, "."}' | cat "${this_hdr}" - | sed 's/^0chrom/chrom/' > "${outprefix_cleanIntronDels}" || true
else
  echo 'Use the temporary file of only clean_intron_deletions as a reference for this sample.'
  echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${tmpfile_cleanIntronDels}" '| grep -v ^chrom | cat' "${this_hdr}" '- | cut -d$\t -f"1-$cut1,$cut1_plus_4" | \'
  echo '  sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_cleanIntronDels}"
  bedtools intersect -wao -a "${this_input}" -b "${tmpfile_cleanIntronDels}" | grep -v '^chrom' | cat "${this_hdr}" - | cut -d$'\t' -f"1-$cut1,$cut1_plus_4" | \
    sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_cleanIntronDels}" || true
  echo ''
fi

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_cleanIntronDels}" '>' "${outprefix_cleanIntronDels_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_cleanIntronDels}" > "${outprefix_cleanIntronDels_collapsed}"
echo ''
this_input="${outprefix_cleanIntronDels_collapsed}"



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with refseq introns and identify refseq introns that are now too short for the spliceosome.'
echo ''
echo 'Bryen SJ, Joshi H, Evesson FJ, Girard C, Ghaoui R, Waddell LB, … Cooper ST (2019) Pathogenic Abmornal Splicing due to Intronic Deletions that induce Biophysical Space Constraint for Spliceosome Assembly. The American Journal of Human Genetics. 105(3), 573-587. Doi: 10.1016/j.ajhg.2019.07.103.'
echo "'We therefore advocate scrutiny of any deletion in a phenotypically consistent gene that renders overall intron length <71 nt (0.1th percentile among human introns) or 5'SS-branchpoint length <59 nt (0.1th percentile among human introns), with our data alerting extreme risk for splicing abnormalities for introns with 5'SS-branchpoint length reduced to <50 nt.'"
echo ''

# head 
# chr1    67093605        67096250        C1orf141        -
# chr1    67096322        67103236        C1orf141        -
# chr1    33540699        33541128        CSMD2   -
# chr1    33541310        33542718        CSMD2   -

num_input_cols=$(head -n 1 $this_input | sed -e 's/\t/\n/g' | wc -l)
num_ref_cols=$(head -n 1 $ucsc_refseq_introns | sed -e 's/\t/\n/g' | wc -l)
overlap_col=$(( num_input_cols + num_ref_cols + 1 ))
overlap_col_plus_1=$(( overlap_col + 1 ))
ref_start_col=$(( num_input_cols + 2 ))
ref_end_col=$(( num_input_cols + 3 ))
svtype_col=-1
cols_string=$(head -n 1 $this_input)
IFS=$'\t' read -r -a array <<< "$cols_string"
i=0
for this_col in "${array[@]}"; do
  i=$(( i + 1 ))
  if [[ $this_col == "SVTYPE" ]]; then
    svtype_col=$i
  fi
  if [[ $this_col == "cnv_type" ]]; then
    svtype_col=$i
  fi
done

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".cleanIntronDels_header."${random_number}".txt
hdr_bit1=`head -n 1 "${this_input}"`
echo -e "${hdr_bit1}\tintron_too_short" > "${this_hdr}"

echo 'bedtools intersect -a' $this_input '-b' $ucsc_refseq_introns '-wao | \'
echo '  awk -v svtype_col='"$svtype_col" '-v overlap_col='"$overlap_col" '-v ref_start_col='"$ref_start_col" '-v ref_end_col='"$ref_end_col" 'BEGIN {FS="\t";OFS="\t"} {if ((($svtype_col=="DEL") || ($svtype_col=="INDEL") || ($svtype_col=="deletion")) && ($ref_start_col<=$2) && ($ref_end_col>=$3) && ((($ref_end_col-$ref_start_col)-($3-$2)) < 70)) {result_intron = ($ref_end_col-$ref_start_col)-($3-$2); print $0, "intron_too_short_for_spliceosome__" result_intron} else {print $0, "."}}' '| \'
echo '  cut -d$''\t' '-f"1-'$num_input_cols','$overlap_col_plus_1'" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- >' "${outprefix_intronTooShort}" '|| true'

bedtools intersect -a $this_input -b $ucsc_refseq_introns -wao | \
  awk -v svtype_col="$svtype_col" -v overlap_col="$overlap_col" -v ref_start_col="$ref_start_col" -v ref_end_col="$ref_end_col" 'BEGIN {FS="\t";OFS="\t"} {if ((($svtype_col=="DEL") || ($svtype_col=="INDEL") || ($svtype_col=="deletion")) && ($ref_start_col<=$2) && ($ref_end_col>=$3) && ((($ref_end_col-$ref_start_col)-($3-$2)) < 70)) {result_intron = ($ref_end_col-$ref_start_col)-($3-$2); print $0, "intron_too_short_for_spliceosome__" result_intron} else {print $0, "."}}' | \
  cut -d$'\t' -f"1-$num_input_cols,$overlap_col_plus_1" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat "${this_hdr}" - > "${outprefix_intronTooShort}" || true

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row.awk' "${outprefix_intronTooShort}" '>' "${outprefix_intronTooShort_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_intronTooShort}" > "${outprefix_intronTooShort_collapsed}"
echo ''
this_input="${outprefix_intronTooShort_collapsed}"



if [[ "$ucsc_gencode_introns" != "" ]]; then

echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with gencode introns and identify gencode introns that are now too short for the spliceosome.'
echo ''
echo 'Bryen SJ, Joshi H, Evesson FJ, Girard C, Ghaoui R, Waddell LB, ... Cooper ST (2019) Pathogenic Abmornal Splicing due to Intronic Deletions that induce Biophysical Space Constraint for Spliceosome Assembly. The American Journal of Human Genetics. 105(3), 573-587. Doi: 10.1016/j.ajhg.2019.07.103.'
echo "'We therefore advocate scrutiny of any deletion in a phenotypically consistent gene that renders overall intron length <71 nt (0.1th percentile among human introns) or 5'SS-branchpoint length <59 nt (0.1th percentile among human introns), with our data alerting extreme risk for splicing abnormalities for introns with 5'SS-branchpoint length reduced to <50 nt.'"
echo ''

# head 
# chr1    67093605        67096250        C1orf141        -
# chr1    67096322        67103236        C1orf141        -
# chr1    33540699        33541128        CSMD2   -
# chr1    33541310        33542718        CSMD2   -

num_input_cols=$(head -n 1 $this_input | sed -e 's/\t/\n/g' | wc -l)
num_ref_cols=$(head -n 1 $ucsc_gencode_introns | sed -e 's/\t/\n/g' | wc -l)
overlap_col=$(( num_input_cols + num_ref_cols + 1 ))
overlap_col_plus_1=$(( overlap_col + 1 ))
ref_start_col=$(( num_input_cols + 2 ))
ref_end_col=$(( num_input_cols + 3 ))
svtype_col=-1
cols_string=$(head -n 1 $this_input)
IFS=$'\t' read -r -a array <<< "$cols_string"
i=0
for this_col in "${array[@]}"; do
  i=$(( i + 1 ))
  if [[ $this_col == "SVTYPE" ]]; then
    svtype_col=$i
  fi
  if [[ $this_col == "cnv_type" ]]; then
    svtype_col=$i
  fi
done

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".cleanIntronDels_header."${random_number}".txt
hdr_bit1=`head -n 1 "${this_input}"`
echo -e "${hdr_bit1}\tgencode_intron_too_short" > "${this_hdr}"

echo 'bedtools intersect -a' $this_input '-b' $ucsc_gencode_introns '-wao | \'
echo '  awk -v svtype_col='"$svtype_col" '-v overlap_col='"$overlap_col" '-v ref_start_col='"$ref_start_col" '-v ref_end_col='"$ref_end_col" 'BEGIN {FS="\t";OFS="\t"} {if ((($svtype_col=="DEL") || ($svtype_col=="INDEL") || ($svtype_col=="deletion")) && ($ref_start_col<=$2) && ($ref_end_col>=$3) && ((($ref_end_col-$ref_start_col)-($3-$2)) < 70)) {result_intron = ($ref_end_col-$ref_start_col)-($3-$2); print $0, "intron_too_short_for_spliceosome__" result_intron} else {print $0, "."}}' '| \'
echo '  cut -d$''\t' '-f"1-'$num_input_cols','$overlap_col_plus_1'" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat' "${this_hdr}" '- >' "${outprefix_gencode_intronTooShort}" '|| true'

bedtools intersect -a $this_input -b $ucsc_gencode_introns -wao | \
  awk -v svtype_col="$svtype_col" -v overlap_col="$overlap_col" -v ref_start_col="$ref_start_col" -v ref_end_col="$ref_end_col" 'BEGIN {FS="\t";OFS="\t"} {if ((($svtype_col=="DEL") || ($svtype_col=="INDEL") || ($svtype_col=="deletion")) && ($ref_start_col<=$2) && ($ref_end_col>=$3) && ((($ref_end_col-$ref_start_col)-($3-$2)) < 70)) {result_intron = ($ref_end_col-$ref_start_col)-($3-$2); print $0, "intron_too_short_for_spliceosome__" result_intron} else {print $0, "."}}' | \
  cut -d$'\t' -f"1-$num_input_cols,$overlap_col_plus_1" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | cat "${this_hdr}" - > "${outprefix_gencode_intronTooShort}" || true

# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row.awk' "${outprefix_gencode_intronTooShort}" '>' "${outprefix_gencode_intronTooShort_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_gencode_intronTooShort}" > "${outprefix_gencode_intronTooShort_collapsed}"
echo ''
this_input="${outprefix_gencode_intronTooShort_collapsed}"

fi



echo ''
echo 'strucVars_annotate_multiple_bedtools_hits.pbs: bedtools intersect with hotspot exons'
echo ''

# #Chr    Start   End     hotspot_exons
# chr1    939274  939412  -1
# chr1    939274  939460  -1
# chr1    1050232 1050329 0.636561334
# chr1    1312786 1312949 0.533326328
# chr1    1564772 1564916 0.715341866

this_hdr="${tmpdir}"/temp_hdr."${sample}"."${script_name}".hotspot_exons."${random_number}".txt

cut1=$(head -n 1 "${this_input}" | sed -e 's/\t/\n/g' | wc -l)
cut1_plus_4=$((cut1+4))
echo 'head -n 1' "${this_input}" '| sed s/^chrom/0chrom/ | sed -e s/$/\t.\t.\t.\thotspot_exons/ >' "${this_hdr}"
head -n 1 "${this_input}" | sed 's/^chrom/0chrom/' | sed -e 's/$/\t.\t.\t.\thotspot_exons/' > "${this_hdr}"

echo 'bedtools intersect -wao -a' "${this_input}" '-b' "${hotspot_exons}" '| grep -v ^chrom | cat' "${this_hdr}" '- | \'
echo '  cut -d$\t -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed s/^0chrom/chrom/ >' "${outprefix_hotspot_exons}"
bedtools intersect -wao -a "${this_input}" -b "${hotspot_exons}" | grep -v '^chrom' | cat "${this_hdr}" - | \
  cut -d$'\t' -f"1-$cut1,$cut1_plus_4" | sort | uniq | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | sed 's/^0chrom/chrom/' > "${outprefix_hotspot_exons}" || true
rm "${this_hdr}"
echo ''

num_hotspot_cols=$(head -n 1 "${hotspot_exons}" | cut -d$'\t' -f4- | sed -e 's/\t/\n/g' | wc -l)
num_output_cols=$(head -n 1 "${outprefix_hotspot_exons}" | sed -e 's/\t/\n/g' | wc -l)
list_columns_to_collapse=""
i=$(( $num_output_cols - $num_hotspot_cols + 1 ))
while [[ $i -le $num_output_cols ]]; do
  if [[ $list_columns_to_collapse == "" ]]; then
    list_columns_to_collapse="${i}"
  else
    new_list_columns_to_collapse="${list_columns_to_collapse},${i}"
    list_columns_to_collapse=$new_list_columns_to_collapse
  fi
  i=$(( $i + 1 ))
done
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column
echo 'awk -f' $sw'/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk' "${outprefix_hotspot_exons}" '>' "${outprefix_hotspot_exons_collapsed}"
awk -f $sw/victorchang_scripts/collapse_manta_bedtools_SVs_to_unique_row_and_value.awk "${outprefix_hotspot_exons}" > "${outprefix_hotspot_exons_collapsed}"
echo ''
this_input="${outprefix_hotspot_exons_collapsed}"



echo ''
echo 'Output:' $this_input
echo ''



rm $outprefix_sorted
rm $outprefix_CDSexons
rm $outprefix_CDSexons_collapsed
rm $outprefix_genes
rm $outprefix_genes_collapsed
rm $outprefix_bndInsideExon
if [[ "$ucsc_gencode_cdsexons" != "" ]]; then
  rm $outprefix_gencode_CDSexons
  rm $outprefix_gencode_CDSexons_collapsed
  rm $outprefix_gencode_genes
  rm $outprefix_gencode_genes_collapsed
  rm $outprefix_gencode_bndInsideExon
fi
rm $outprefix_EHRFr104
rm $outprefix_EHRFr104_collapsed
rm $outprefix_gnomadsv
rm $outprefix_gnomadsv_AF
rm $outprefix_gnomadsv_collapsed
rm $outprefix_dgv
rm $outprefix_dgv_collapsed
rm $outprefix_segdup
rm $outprefix_dbscSNV_largeSV
rm $outprefix_dbscSNV
rm $outprefix_dbscSNV_1lastCol
rm $outprefix_dbscSNV_collapsed
rm $outprefix_FANTOM5_TSS
rm $outprefix_FANTOM5_TSS_collapsed
rm $outprefix_pext_artery_aorta
rm $outprefix_pext_artery_aorta_collapsed
rm $outprefix_pext_artery_coronary
rm $outprefix_pext_artery_coronary_collapsed
rm $outprefix_pext_heart_atrial_appendage
rm $outprefix_pext_heart_atrial_appendage_collapsed
rm $outprefix_pext_heart_left_ventricle
rm $outprefix_pext_heart_left_ventricle_collapsed
rm $outprefix_pext_genes
rm $outprefix_pext_genes_collapsed
rm $outprefix_cleanIntronDels
rm $outprefix_intronTooShort
if [[ "$ucsc_gencode_introns" != "" ]]; then
  rm $outprefix_intronTooShort_collapsed
  rm $outprefix_gencode_intronTooShort
  rm $outprefix_gencode_intronTooShort_collapsed
else
  rm $outprefix_intronTooShort_collapsed
fi
rm $outprefix_hotspot_exons
# The output file is $outprefix_hotspot_exons_collapsed



echo ''
echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"


