#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=2GB
#PBS -l jobfs=2GB
#PBS -l ncpus=1
#PBS -N verifybamid
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

module load htslib/1.9
module load samtools/1.10

. "${sw_and_refs}"

/g/data/jb96/software/verifybamid2/VerifyBamID/bin/VerifyBamID \
  --SVDPrefix $verifybamid_population_frequencies \
  --Reference "${ref_fasta_fa}" \
  --BamFile "${infile}" > "${outfile}"

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''

# Example output from https://github.com/Griffan/VerifyBamID
#
# Estimation from OptimizeHeter:
# Contaminating Sample PC1:-0.623602      PC2:0.57292
# Intended Sample  PC1:-0.036304  PC2:0.0200112
# Alpha:0.0013662
#
# First line: Which model used
# Second line: PC coordinates of Contaminating Sample
# Third line: PC coordinates of Intended Sample(the sample your are interested)
# Fourth line: Estimated contamination level
