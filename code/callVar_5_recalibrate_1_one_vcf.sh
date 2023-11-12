#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=03:00:00
#PBS -l mem=12G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N VariantRecalibrator
#PBS -lstorage=gdata/abcd

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

module load python3-as-python
module load java/jdk-8.40
module load R/3.6.1

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# indel annotation & tranch values obtained from:
# https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.hg38.wgs.inputs.json
indel_an_values=( FS ReadPosRankSum MQRankSum QD SOR DP )
# convert to formatted string for gatk command
indel_an_values_str=$( for an in ${indel_an_values[@]} ; do printf " -an $an" ; done )
indel_tranche_values=( 100.0 99.95 99.9 99.5 99.0 97.0 96.0 95.0 94.0 93.5 93.0 92.0 91.0 90.0 )
# convert to formatted string for gatk command
indel_tranche_values_str=$( for t in ${indel_tranche_values[@]} ; do printf " -tranche $t" ; done )

# SNP annotation & tranch values obtained from:
# https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/joint-discovery-gatk4-local.hg38.wgs.inputs.json
snp_tranche_values=( 100.0 99.95 99.9 99.8 99.6 99.5 99.4 99.3 99.0 98.0 97.0 90.0 )
snp_tranche_values_str=$( for t in ${snp_tranche_values[@]} ; do printf " -tranche $t" ; done )
snp_an_values=( QD MQRankSum ReadPosRankSum FS MQ SOR DP )
snp_an_values_str=$( for an in ${snp_an_values[@]} ; do printf " -an $an" ; done )

# NOTE that the "--resource" format has changed between versions
# from GATK 4.1.1.0 the format changes to:
# "--resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf}"
# previously was: (note the ":" swaps with the " ")
# "--resource mills,known=false,training=true,truth=true,prior=12:${mills_resource_vcf}"

echo ''
echo 'Run gatk indel VariantRecalibrator' $infile 'to produce output file' $outfile_indel_recal
echo ''

axiomPoly_string="--resource:axiomPoly,known=false,training=true,truth=false,prior=10 "${axiomPoly_resource_vcf}
if [[ $genome_version == "hg19" ]]; then
  axiomPoly_string=""
fi

${gatk_path} --java-options "-server -Xms1g -Xmx8g -Djava.io.tmpdir=$tmpdir" \
    VariantRecalibrator \
    -V ${infile} \
    -O ${outfile_indel_recal} \
    --tranches-file ${outfile_indel_tranches} \
    --rscript-file ${outfile_indel_rscript} \
    --trust-all-polymorphic \
    ${indel_tranche_values_str} \
    ${indel_an_values_str} \
    -mode INDEL \
    --max-gaussians 4 \
    --resource:mills,known=false,training=true,truth=true,prior=12 ${mills_resource_vcf} \
    $axiomPoly_string \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp_resource_vcf}

echo ''
echo 'Run gatk SNP VariantRecalibrator' $infile 'to produce output file' $outfile_snp_recal
echo ''
# Most samples used     --max-gaussians 6 \
# Larger samples might need    --max-gaussians 4 \

${gatk_path} --java-options "-server -Xms1g -Xmx8g -Djava.io.tmpdir=$tmpdir" \
    VariantRecalibrator \
    -V ${infile} \
    -O ${outfile_snp_recal} \
    --tranches-file ${outfile_snp_tranches} \
    --trust-all-polymorphic \
    --rscript-file ${outfile_snp_rscript} \
    ${snp_tranche_values_str} \
    ${snp_an_values_str} \
    -mode SNP \
    --max-gaussians 6 \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap_resource_vcf} \
    --resource:omni,known=false,training=true,truth=true,prior=12 ${omni_resource_vcf} \
    --resource:1000G,known=false,training=true,truth=false,prior=10 ${one_thousand_genomes_resource_vcf} \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp_resource_vcf}

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

