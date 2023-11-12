#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l storage=gdata/abcd
#PBS -l jobFs=32GB
#NOTE - larger inputs (like genome v3.1) require jobfs for sorting
set -o pipefail
# This file converts the preprocessed gnomad VCF file to 
# annovar db format. Processes by chromosome.
#
# Longest job for gnomad v3.1 genomes < 25min
#
: ${chr:? missing variable - chr}
module load bcftools
module load htslib
exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "TERMINATING.  $1" >&2
        exit 1
    fi
}
check_file_exists() {
    if ! [ -f "$1" ] ; then
        echo "TERMINATING.  Unable to find file: $1" >&2
        exit 1
    fi
}
scripts_d=/my/directory/variant_calling_and_annotation_pipeline/annovar
d=/my/directory/variant_calling_and_annotation_pipeline/annovar/humandb
cd "$d"
exit_on_error "Unable to change to dir: $d"

# set run version and other identifying data
exomes_211=exomes_211
genomes_311=genomes_311
genomes_31=genomes_311
hg38_exomes_211=hg38_$exomes_211
hg19_exomes_211=hg19_$exomes_211
hg38_genomes_311=hg38_$genomes_311
hg38_genomes_31=hg38_$genomes_31
#TODO - set run_ver to appropriate for your run
#run_ver=$hg19_exomes_211   # completed
#run_ver=$hg38_exomes_211
run_ver=$hg38_genomes_31

if [ "$run_ver" == "$hg38_exomes_211" ] ; then
    in_file=gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38_processed.vcf
    db_pre=gnomad211_exome_
elif [ "$run_ver" == "$hg19_exomes_211" ] ; then
    in_file=gnomad.exomes.r2.1.1.sites.${chr}.processed.vcf
    db_pre=gnomad211_exome_
elif [ "$run_ver" == "$hg38_genomes_31" ] ; then
    in_file=gnomad.genomes.v3.1.sites.chr${chr}.processed.vcf
    db_pre=gnomad31_genome_
elif [ "$run_ver" == "$hg38_genomes_311" ] ; then
    in_file=gnomad.genomes.v3.1.1.sites.chr${chr}.processed.vcf
    db_pre=gnomad311_genome_
fi
################
check_file_exists "$in_file"
vt_decomp=${in_file}"_vt_decomp.vcf" # 
vt_decomp_blocksub=${in_file}"_vt_decomp_blocksub.vcf" # 
vt_decomp_blocksub_norm=${in_file}"_vt_decomp_blocksub_norm.vcf" # 
final_output_annovar_db=$in_file"_annovar.txt" # 
tmp_output_annovar_db=tmp_$in_file"_annovar.txt" # 
#
#
# TODO: change ref if working with GRCh37
if [ "${run_ver%%_*}" == "hg38" ] ; then
    ref_fasta=/my/directory/variant_calling_and_annotation_pipeline/reference_data/hg38.noalt.decoy.bwa/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
else
    ref_fasta=/my/directory/variant_calling_and_annotation_pipeline/reference_data/hs37d5/hs37d5x.fa
fi
echo "Using ref_fasta: $ref_fasta" >> /dev/stderr
check_file_exists "$ref_fasta"
sw=/my/directory/variant_calling_and_annotation_pipeline
vt=$sw/vt/vt
check_file_exists "$vt"
$vt decompose -s "$in_file" -o "$vt_decomp" ##
exit_on_error "Problem with vt decompose"
bgzip -f "${vt_decomp}"
exit_on_error "problem with bgzip for: $vt_decomp"
tabix -f "${vt_decomp}.gz"
exit_on_error "Poblem with tabix for ${vt_decomp}.gz"
$vt decompose_blocksub -a "${vt_decomp}.gz" -o "$vt_decomp_blocksub"  #
exit_on_error "Problem with vt decompose_blocksub"
bgzip -f "${vt_decomp_blocksub}"
exit_on_error "problem with bgzip for: $vt_decomp_blocksub"
tabix -f "${vt_decomp_blocksub}.gz"
exit_on_error "Poblem with tabix for ${vt_decomp_blocksub}.gz"
$vt normalize -r "$ref_fasta" -o "$vt_decomp_blocksub_norm" "${vt_decomp_blocksub}.gz" 
exit_on_error "vt Normalize error."
#tabix -p vcf "$vt_decomp_blocksub_norm"
#
#
convert2annovar=$sw/annovar/convert2annovar_ed.pl
check_file_exists "$convert2annovar"
awk_script=$scripts_d/extract_final_db_fields.awk
check_file_exists "$awk_script"
perl "$convert2annovar" -format vcf4 "$vt_decomp_blocksub_norm" -includeinfo | \
    awk -v db_field_prefix="$db_pre" -v script_name="$awk_script" \
    -f "$awk_script"  > "$tmp_output_annovar_db"
# sort output & remove duplicate entries
{
    head -1 "$tmp_output_annovar_db"
    tail -n+2 "$tmp_output_annovar_db" | LC_ALL=C sort -V -u 
} > "$final_output_annovar_db"
exit_on_error "Problem with convert2annovar"
#
#
# tidy up
rm "${vt_decomp}.gz" "${vt_decomp}.gz.tbi" "${vt_decomp_blocksub}.gz" "${vt_decomp_blocksub}.gz.tbi" "$vt_decomp_blocksub_norm"
rm "$tmp_output_annovar_db"
echo "SUCCESS" >>/dev/stderr
