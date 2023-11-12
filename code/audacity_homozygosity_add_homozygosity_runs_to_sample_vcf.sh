#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=100G
#PBS -l ncpus=1
#PBS -l jobfs=100G
#PBS -N homannovar
#PBS -lstorage=gdata/abcd

# mem=50G is enough for 900000 vars. more mem needed for 1000000 vars.

set -euo pipefail

#sample=$1
#infile_vcf=$2 # vcf file to be annotated
#infile2_flag=$3 # BigHomozygRun
#infile2=$4 # bed file containing BigHomozygRun runs of homozygosity
#infile3_flag=$5 # HRun
#infile3=$6 # bed file containing HRun runs of homozygosity
#outdir=$7
#outfile_basename=$8
#sw_and_refs=$9
#sample_job_differentiating_id="${10}"

outfile="${outdir}"/"${outfile_basename}"
outfile_prefix="${outfile%.gz}"
outfile_prefix="${outfile_prefix%.vcf}"

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

module load htslib
module load bcftools
module load python3
module load R/3.6.1
module load intel-compiler/2019.3.199


. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

##################################################
# Convert the audacity homozygosity result to an annovar regions annotation table
##################################################

tmp_annovar2="${tmpdir}"/hg38_large_audacity_homozygosity_runs.txt
tmp_annovar3="${tmpdir}"/hg38_audacity_homozygosity_runs.txt

echo -e $infile2_flag"\tchr1\t1\t2\t"$infile2_flag > $tmp_annovar2 # make sure got at least 1 record so annovar wont crash
echo 'awk -v flag='"$infile2_flag" 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print flag, $1, $2, $3, flag}}' $infile2 '| sort -k1,1V -k2,2V -k3,3V | uniq >>' $tmp_annovar2
awk -v flag="$infile2_flag" 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print flag, $1, $2, $3, flag}}' $infile2 | sort -k1,1V -k2,2V -k3,3V | uniq >> $tmp_annovar2
echo ''

echo -e $infile3_flag"\tchr1\t1\t2\t"$infile3_flag > $tmp_annovar3 # make sure got at least 1 record so annovar wont crash
echo 'awk -v flag='"$infile3_flag" 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print flag, $1, $2, $3, flag}}' $infile3 '| sort -k1,1V -k2,2V -k3,3V | uniq >>' $tmp_annovar3
awk -v flag="$infile3_flag" 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print flag, $1, $2, $3, flag}}' $infile3 | sort -k1,1V -k2,2V -k3,3V | uniq >> $tmp_annovar3
echo ''

##################################################
# Annotate the SNP/indel VCF with the audacity homozygosity regions
##################################################

outfile_basename=$(basename $outfile)
tmp_annotated="${tmpdir}"/"${outfile_basename}"

echo $sw'/annovar/table_annovar.pl' "$infile_vcf" ${tmpdir}'/ -vcfinput -buildver' $genome_version '\'
echo '  -out' "$outfile_prefix" '-remove \'
echo '  -protocol large_audacity_homozygosity_runs,audacity_homozygosity_runs \'
echo '  -operation r,r -nastring . \'
echo '  -arg -time,-time'

$sw/annovar/table_annovar.pl "$infile_vcf" ${tmpdir}/ -vcfinput -buildver "$genome_version" \
  -out "$outfile_prefix" -remove \
  -protocol large_audacity_homozygosity_runs,audacity_homozygosity_runs \
  -operation r,r -nastring . \
  -arg -time,-time

echo ''

echo 'mv' "${outfile_prefix}"'."${genome_version}"_multianno.vcf' "${outfile_prefix}"'.vcf'
mv "${outfile_prefix}"."${genome_version}"_multianno.vcf "${outfile_prefix}".vcf
echo ''

echo ''
echo 'outfile:' "${outfile_prefix}".vcf
echo ''

echo 'Finished!'
echo ''

touch "${done_file}"
rm -f "${lock_file}"



