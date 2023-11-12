#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=20:00:00
#PBS -l mem=64G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N vep
#PBS -lstorage=gdata/abcd

#set -euo pipefail

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
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

PATH="/g/data/jb96/software/vep/local_lib/perl5/bin${PATH:+:${PATH}}"; export PATH;
PERL5LIB="/g/data/jb96/software/vep/local_lib/perl5/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB
PERL_LOCAL_LIB_ROOT="/g/data/jb96/software/vep/local_lib/perl5${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"; export PERL_LOCAL_LIB_ROOT;
PERL_MB_OPT="--install_base \"/g/data/jb96/software/vep/local_lib/perl5\""; export PERL_MB_OPT;
PERL_MM_OPT="INSTALL_BASE=/g/data/jb96/software/vep/local_lib/perl5"; export PERL_MM_OPT;

export PERLBREW_ROOT=/g/data/jb96/software/vep/perlbrew_installation
source /g/data/jb96/software/vep/perlbrew_installation/etc/bashrc
perlbrew use perl-5.31.10
module load htslib
module load samtools/1.10
module load python3
export PERL5LIB=$PERL5LIB:/g/data/jb96/software/vep/vep_cache/Plugins
export PERL5LIB=$PERL5LIB:$vep_loftee
export PERL5LIB=$PERL5LIB:/g/data/jb96/software/vep_plugins/UTRannotator/UTRannotator

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

outfile2="${outfile%.gz}"
outfile2="${outfile2%.vcf}"
outfile2="${outfile2%.txt}"
outfile2="${outfile2%.vep}".vep_reformatted.vcf
rm -rf $outfile2

# VEP: Annotate with basic VEP, not everything.
# VEP: Annotate with VEP plugins loftee and UTRannotator, not with all the VEP plugins.

echo 'perl /g/data/jb96/software/vep/ensembl-vep/vep --cache --merged --dir_cache /g/data/jb96/software/vep/vep_cache --dir_plugins /g/data/jb96/software/vep/ensembl-vep \'
echo '  --offline -i' $infile '-o' $outfile '--use_given_ref --force_overwrite --vcf --format vcf --assembly' $vep_assembly '--regulatory --canonical --hgvs --hgvsg --fasta' $ref_fasta_fa '\'
echo '  --plugin TSSDistance \'
echo '  --plugin miRNA \'
echo '  --plugin LoFtool,/g/data/jb96/software/vep_plugins/VEP_plugins/VEP_plugins/LoFtool_scores.txt \'
echo '  --plugin MPC,/g/data/jb96/software/vep_plugins/MPC/fordist_constraint_official_mpc_values_v2.txt.gz \'
echo '  --plugin UTRannotator,'$vep_utrannotator '\'
echo '  --plugin LoF,loftee_path:'$vep_loftee',human_ancestor_fa:'$vep_loftee_human_ancestor_fa',conservation_file:'$vep_loftee_conservation_file
echo ''

perl /g/data/jb96/software/vep/ensembl-vep/vep --cache --merged --dir_cache /g/data/jb96/software/vep/vep_cache --dir_plugins /g/data/jb96/software/vep/ensembl-vep \
  --offline -i $infile -o $outfile --use_given_ref --force_overwrite --vcf --format vcf --assembly $vep_assembly --regulatory --canonical --hgvs --hgvsg --fasta $ref_fasta_fa \
  --plugin TSSDistance \
  --plugin miRNA \
  --plugin LoFtool,/g/data/jb96/software/vep_plugins/VEP_plugins/VEP_plugins/LoFtool_scores.txt \
  --plugin MPC,/g/data/jb96/software/vep_plugins/MPC/fordist_constraint_official_mpc_values_v2.txt.gz \
  --plugin UTRannotator,$vep_utrannotator \
  --plugin LoF,loftee_path:$vep_loftee:$vep_loftee_human_ancestor_fa,$vep_loftee_conservation_file

echo ''

set -euo pipefail

echo ''
echo 'python3' $sw'/victorchang_scripts/convert_vep_multiple_transcript_vcf_to_official_vcf_format.py -i' "${outfile}" '-o' "${outfile2}"
echo ''

python3 $sw/victorchang_scripts/convert_vep_multiple_transcript_vcf_to_official_vcf_format.py -i "${outfile}" -o "${outfile2}"

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''
