#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=40G
#PBS -N annovar
#PBS -lstorage=gdata/abcd

set -euo pipefail

#infile=$1 # vcf input
#outdir=$2
#outfile=$3 # vcf output
#chrom=$4
#sw_and_refs=$5

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
if [[ $genome_version == "hg38" ]]; then
  echo 'Get the hg19 coordinates for this hg38 coordinates, using liftover.'
fi
if [[ $genome_version == "hg19" ]]; then
  echo 'Get the hg38 coordinates for this hg19 coordinates, using liftover.'
fi
echo ''

tmp_hg19_in_bed="${tmpdir}"/tmp_in_liftover.hg19_liftover.chrom_"${chrom}".bed
tmp_hg19_out_bed="${tmpdir}"/tmp_out_liftover.hg19_liftover.chrom_"${chrom}".bed
tmp_hg19_rejected_bed="${tmpdir}"/tmp_out_rejected_liftover.hg19_liftover.chrom_"${chrom}".bed
tmp_hg19_out_vcf="${tmpdir}"/tmp_out_vcf.hg19_liftover.chrom_"${chrom}".vcf

# create input - will look like:
#chrY	2781760	2781761	chrY|2781761|CA|C
#chrY	2781760	2781761	chrY|2781761|CAA|C

echo "Creating bed file from VCF for liftover"
echo ''

#echo 'bcftools query -f' '%CHROM\t%POS0\t%POS\t%CHROM|%POS|%REF|%ALT\n' $infile '>' $tmp_hg19_in_bed
#bcftools query -f '%CHROM\t%POS0\t%POS\t%CHROM|%POS|%REF|%ALT\n' $infile '>' $tmp_hg19_in_bed
#echo ''

# The input to liftover must have 'chr' in front of chromosome, even if it is hg19 without chr
echo 'grep -v' '^#' $infile '| awk' 'BEGIN {FS="\t";OFS="\t"} {print $1, $2-1, $2, $1 "|" $2 "|" $4 "|" $5}' '| sed' 's/^/chr/' '| sed' 's/^chrchr/chr/' '>' $tmp_hg19_in_bed
grep -v '^#' $infile | awk 'BEGIN {FS="\t";OFS="\t"} {print $1, $2-1, $2, $1 "|" $2 "|" $4 "|" $5}' | sed 's/^/chr/' | sed 's/^chrchr/chr/' > $tmp_hg19_in_bed
echo ''

# execute the liftover
# Successful liftover file like:
#chrY	2649801	2649802	chrY|2781761|CA|C
#chrY	2649801	2649802	chrY|2781761|CAA|C

echo "Performing liftover"
echo ''

echo $sw'/liftover/liftOver' $tmp_hg19_in_bed $liftover_chain $tmp_hg19_out_bed $tmp_hg19_rejected_bed
$sw/liftover/liftOver $tmp_hg19_in_bed $liftover_chain $tmp_hg19_out_bed $tmp_hg19_rejected_bed
echo ''

if [[ $genome_version == "hg38" ]]; then
  echo 'Adding hg19 annotation to VCF INFO field'
  to_genome=hg19
fi
if [[ $genome_version == "hg19" ]]; then
  echo 'Adding hg38 annotation to VCF INFO field'
  to_genome=hg38
fi

echo 'awk -f' $sw'/victorchang_scripts/add_liftover_results_as_annotation_to_vcf.awk -v to_genome='"$to_genome" '-v liftover_results_bed='"$tmp_hg19_out_bed" $infile '>' $tmp_hg19_out_vcf
awk -f $sw/victorchang_scripts/add_liftover_results_as_annotation_to_vcf.awk -v to_genome="$to_genome" -v liftover_results_bed="$tmp_hg19_out_bed" $infile > $tmp_hg19_out_vcf
echo ''

if [[ $gnomad_constraints == "" ]]; then

  echo 'cp' $tmp_hg19_out_vcf $outfile
  cp $tmp_hg19_out_vcf $outfile
  echo ''

else

  tmp_gnomad_in_vcf_bed="${tmpdir}"/tmp_in_vcf_bed.gnomad_intersect.chrom_"${chrom}".bed
  tmp_gnomad_in_gnomad_bed="${tmpdir}"/tmp_in_gnomad_bed.gnomad_intersect.chrom_"${chrom}".bed
  tmp_gnomad_intersect_bed="${tmpdir}"/tmp_in_intersect_bed.gnomad_intersect.chrom_"${chrom}".bed

  ################################ >>
  # STEP 1 - create input bed from VCF - will look like:
  #chrY	2781760	2781761	chrY|2781761|CA|C
  #chrY	2781760	2781761	chrY|2781761|CAA|C

  echo "Creating bed file from VCF input"
  echo ''

  # NOTE there are potentially duplicate entries after the vt software splits 
  # multi-allelic entries
  # NOTE embedded tab in sort -t option
  # NOTE uniq -f3 to ignore first 3 fields in determining if uniq
  echo 'bcftools query -f' '%CHROM\t%POS\t%END\t%CHROM|%POS|%REF|%ALT\n' $tmp_hg19_out_vcf '| LC_ALL=C sort -t"	" -k4,4 | uniq -f3 >' $tmp_gnomad_in_vcf_bed
  bcftools query -f '%CHROM\t%POS\t%END\t%CHROM|%POS|%REF|%ALT\n' $tmp_hg19_out_vcf | LC_ALL=C sort -t"	" -k4,4 | uniq -f3 > $tmp_gnomad_in_vcf_bed
  echo ''

  ################################ >>
  # STEP 2 - create a cut-down version of the gnomad constraints file
  # to use as a bed file to intersect with VCF bed file.
  # Reason why use cutdown version: there are many columns in the gnomad constraints file
  # and we can save space by not getting them until we merge back into the VCF annotation
  # gnomad_constraints header:
  #chrom	start	end	gnomad_constraints_gene	gnomad_constraints_syn_z	gnomad_constraints_oe_syn	gnomad_constraints_oe_syn_lower	gnomad_constraints_oe_syn_upper	gnomad_constraints_mis_z	gnomad_constraints_oe_mis	gnomad_constraints_oe_mis_lower	gnomad_constraints_oe_mis_upper	gnomad_constraints_pLI	gnomad_constraints_oe_lof	gnomad_constraints_oe_lof_lower	gnomad_constraints_oe_lof_upper
  # id will be <chrom>|<start>|<end>|<gnomad_constraints_gene>
  # output will be: <chrom>\t<start>\t<end>\t<id>
  # First check the header is as expected - NOTE converted tabs to spaces
  expected="chrom start end gnomad_constraints_gene"
  header=$(  head -1 "$gnomad_constraints" | cut -f1-4 | tr '\t' ' ' )
  echo "header=$(  head -1 "$gnomad_constraints" | cut -f1-4 | tr '\t' ' ' )"
  echo ''

  if [ "$expected" != "$header" ] ; then
    echo "TERMINATING.  Issue with gnomad constraints file format."
    exit 1
  fi
  echo "Creating abbreviated bed file from gnomad constraints input"
  echo ''

  echo "awk -v chr="$chrom" 'BEGIN {script}'" $gnomad_constraints '>' $tmp_gnomad_in_gnomad_bed
  awk -v chr="$chrom" 'BEGIN {
FS="\t"
OFS = "\t"
}
NR>1 && $1==chr {
# I have found true (all fields same) duplicate entries in the gnomad constraints file
# Do not include them to make things tidier downstream
if (!($0 in seen)) {
    print $1, $2, $3, $1 "|" $2 "|" $3 "|" $4
    seen[$0]++
}
}' $gnomad_constraints > $tmp_gnomad_in_gnomad_bed
  echo ''

  ################################ >>
  # STEP 3 - use bedtools to get overlap of VCF sites with gnomad sites
  echo "Bedtools intersect with gnomad constraints file"
  echo ''
  #NOTE - only interested in matching entries
  echo 'bedtools intersect -wao -a' "$tmp_gnomad_in_vcf_bed" '-b' "$tmp_gnomad_in_gnomad_bed" '| awk -F"\t"' '$5!="."' '>' $tmp_gnomad_intersect_bed
  bedtools intersect -wao -a "$tmp_gnomad_in_vcf_bed" -b "$tmp_gnomad_in_gnomad_bed" | awk -F"\t" '$5!="."' > $tmp_gnomad_intersect_bed
  echo ''

  ################################ >>
  # STEP 4 - Use awk script to append the annotation to the INFO field
  echo "Adding gnomad gene constraint annotations to VCF INFO field" >> /dev/stderr
  echo ''

  echo 'awk -f' $sw'/victorchang_scripts/add_gnomad_gene_constraints_to_vcf.awk -v chr='"$chrom" '-v bed='"$tmp_gnomad_intersect_bed" '-v gnom='"$gnomad_constraints" $tmp_hg19_out_vcf '>' $outfile
  awk -f $sw/victorchang_scripts/add_gnomad_gene_constraints_to_vcf.awk -v chr="$chrom" -v bed="$tmp_gnomad_intersect_bed" -v gnom="$gnomad_constraints" $tmp_hg19_out_vcf > $outfile
  echo ''

fi

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo 'outfile:' $outfile
echo ''
