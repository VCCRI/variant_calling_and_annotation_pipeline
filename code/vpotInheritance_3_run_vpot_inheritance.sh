#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=80G
#PBS -l ncpus=1
#PBS -l jobfs=100G
#PBS -N vpotInheritance
#PBS -lstorage=gdata/abcd

set -euo pipefail

#family=$1
#inheritance_model=$2
#infile=$3
#pedigree_file=$4
#maf=$5
#outdir=$6
#outfile=$7
#sw_and_refs=$8

echo ''
echo 'family:' $family
echo 'inheritance_model:' $inheritance_model
echo 'infile:' $infile
echo 'pedigree_file:' $pedigree_file
echo 'maf:' $maf
echo 'outdir:' $outdir
echo 'outfile:' $outfile
echo 'sw_and_refs:' $sw_and_refs
echo ''

out_prefix=$outfile

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

module load intel-mkl/2019.3.199
module load python3

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

rand_num=$RANDOM
tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

# We will do the following per family:
# AD = autosomal dominant
# AR = autosomal recessive
# CH = compound hets
# DN = denovo
# MAF = 0.01
# This script will have been called to do one of them.

# Input pedigree file contains the true pedigree and phenotypes
# eg.
# 10077	10077	14407	14406	1	2
# 10077	14406			2	1
# 10077	14407			1	1
# eg.
# 2000055	2000055	2000059	2000054	1	2
# 2000055	2000054			2	1
# 2000055	2000056	2000059	2000054	1	1
# 2000055	2000059			1	1
# eg.
# 13560	14048	14047	14046	2	1
# 13560	14049	14047	14046	2	1
# 13560	14050	14047	14046	1	1

# Pedigree files used by this script:

#
# by-families:
# 10077	10077	14407	14406	1	2
# 10077	14406			2	1
# 10077	14407			1	1
#
# father pedigree for AD father:
# 10077	10077	14407	14406	1	2
# 10077	14406			2	1
# 10077	14407			1	2
#
# mother pedigree for AD father:
# 10077	10077	14407	14406	1	2
# 10077	14406			2	2
# 10077	14407			1	1

# Name the various possible output files
out_prefix_other_models="${out_prefix}"."${inheritance_model}"_
out_prefix_AD_father="${out_prefix}"."${inheritance_model}"_father_
out_prefix_AD_mother="${out_prefix}"."${inheritance_model}"_mother_

# Name the various input pedigree files
tmp_pedigree="${tmpdir}"/"${family}".tmp_pedigree.txt
tmp_vpot_AD_father_ped="${tmpdir}"/"${family}".tmp_vpot_AD_father_ped.txt
tmp_vpot_AD_mother_ped="${tmpdir}"/"${family}".tmp_vpot_AD_mother_ped.txt

# Start creating pedigree input files for just this family
grep '^'$family $pedigree_file > $tmp_pedigree

# Work out where the father and mother are, for autosomal dominant processing
trio=$(head -n 1 $tmp_pedigree | sed -e 's/\t/ /g')
IFS=' ' read -r -a array <<< "$trio"
proband="${array[1]}"
father="${array[2]}"
mother="${array[3]}"

# Create a temporary input file that contains just the family members, so that the output will contain just the family members

infile_basename=$(basename $infile)
tmp_hdr="${tmpdir}"/tmp_hdr."${infile_basename}"
tmp_infile1="${tmpdir}"/tmp_infile1."${infile_basename}"
tmp_infile2="${tmpdir}"/tmp_infile2."${infile_basename}"

list_family_members=$(grep -v '^Family' $pedigree_file | grep $family | cut -d$'\t' -f2)

list_family_member_columns=''
for this_sample in $list_family_members; do
  if [[ $this_sample != '0' ]]; then
    this_sample_col=$(head -n 1 $infile | sed -e 's/\t/\n/g' | grep -n $this_sample | cut -d':' -f1 || true)
    if [[ $this_sample_col != '' ]]; then
      if [[ $list_family_member_columns == '' ]]; then
        list_family_member_columns=$this_sample_col
      else
        list_family_member_columns=$list_family_member_columns','$this_sample_col
      fi
    fi
  fi
done

echo 'Extract just the samples in this family, as input to vpot samplef for inheritance.'
echo 'awk -v list_cols='"$list_family_member_columns" 'BEGIN {FS="\t";OFS="\t"} {printf $1; for (i=2; i<=11; i++) {printf OFS $i}; num=split(list_cols, arr, ","); for (i=1; i<=num; i++) {j=arr[i]; printf OFS $j}; printf "\n"}' $infile '>' $tmp_infile1
awk -v list_cols="$list_family_member_columns" 'BEGIN {FS="\t";OFS="\t"} {printf $1; for (i=2; i<=11; i++) {printf OFS $i}; num=split(list_cols, arr, ","); for (i=1; i<=num; i++) {j=arr[i]; printf OFS $j}; printf "\n"}' $infile > $tmp_infile1
echo ''

echo 'Also sort because chroms are in order 1..21 X 22 M Y, and should be in order 1..22 M X Y'
# Input file is VPOL output from vpot priority: Ranking   Priority_Score   #CHROM   POS   ID   REF   ALT
head -n 1 $tmp_infile1 > $tmp_hdr
echo 'awk' 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $tmp_infile1 '| sort -k1,1V -k2,2V -k3,3V -k4,4V -k6,6 -k7,7 | cat' $tmp_hdr '- >' $tmp_infile2
awk 'BEGIN {FS="\t";OFS="\t"} {if (NR>1) {print $0}}' $tmp_infile1 | sort -k1,1V -k2,2V -k3,3V -k4,4V -k6,6 -k7,7 | cat $tmp_hdr - > $tmp_infile2
echo ''

# Now run the vpot inheritance model

if [[ "$inheritance_model" == "AD" ]] ; then

    # For autosomal dominant pedigree input, set one of the parents to having the disease phenotype and the other to unaffected
    awk -v father="$father" -v mother="$mother" 'BEGIN {FS="\t";OFS="\t"} {if ($2 == father) {$6 = "2"}; if ($2 == mother) {$6 = "1"}; print $0}' $tmp_pedigree > $tmp_vpot_AD_father_ped
    awk -v father="$father" -v mother="$mother" 'BEGIN {FS="\t";OFS="\t"} {if ($2 == mother) {$6 = "2"}; if ($2 == father) {$6 = "1"}; print $0}' $tmp_pedigree > $tmp_vpot_AD_mother_ped

    if [[ $father != "ND" ]]; then
        # AD autosomal dominant: father

        echo 'Pedigree used for AD-father:'
        cat $tmp_vpot_AD_father_ped
        echo ''

        echo 'python3' $sw'/VPOT/VPOT.py samplef' $out_prefix_AD_father $tmp_infile2 $tmp_vpot_AD_father_ped $proband $inheritance_model
        python3 $sw/VPOT/VPOT.py samplef $out_prefix_AD_father $tmp_infile2 $tmp_vpot_AD_father_ped $proband $inheritance_model
        echo ''
    fi

    if [[ $mother != "ND" ]]; then
        # AD autosomal dominant: mother

        echo 'Pedigree used for AD-mother:'
        cat $tmp_vpot_AD_mother_ped
        echo ''

        echo 'python3' $sw'/VPOT/VPOT.py samplef' $out_prefix_AD_mother $tmp_infile2 $tmp_vpot_AD_mother_ped $proband $inheritance_model
        python3 $sw/VPOT/VPOT.py samplef $out_prefix_AD_mother $tmp_infile2 $tmp_vpot_AD_mother_ped $proband $inheritance_model
        echo ''
    fi

    if [[ $father != "ND" ]]; then
        echo 'Rename the father vpot output file to remove the random number in the file name.'
        vpot_output_file_pattern="${out_prefix_AD_father}${proband}"_"${inheritance_model}"_variant_filtered_output_file_*.txt
        vpot_output_file_renamed1="${out_prefix_AD_father}${proband}"_"${inheritance_model}"_variant_filtered_output_file.txt
        vpot_output_file=$(ls -t $vpot_output_file_pattern | head -n 1)
        echo 'sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7' $vpot_output_file '>' $vpot_output_file_renamed1
        sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7 $vpot_output_file > $vpot_output_file_renamed1
        echo ''
        rm $vpot_output_file

        echo 'output file:' $vpot_output_file_renamed1
        echo ''
    fi

    if [[ $mother != "ND" ]]; then
        echo 'Rename the mother vpot output file to remove the random number in the file name.'
        vpot_output_file_pattern="${out_prefix_AD_mother}${proband}"_"${inheritance_model}"_variant_filtered_output_file_*.txt
        vpot_output_file_renamed2="${out_prefix_AD_mother}${proband}"_"${inheritance_model}"_variant_filtered_output_file.txt
        vpot_output_file=$(ls -t $vpot_output_file_pattern | head -n 1)
        echo 'sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7' $vpot_output_file '>' $vpot_output_file_renamed2
        sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7 $vpot_output_file > $vpot_output_file_renamed2
        echo ''
        rm $vpot_output_file

        echo 'output file:' $vpot_output_file_renamed2
        echo ''
    fi

elif [[ "$inheritance_model" == "CaseControl" ]] ; then

    cat $tmp_pedigree
    echo ''

    echo 'python3' $sw'/VPOT/VPOT.py samplef' $out_prefix_other_models $tmp_infile2 $tmp_pedigree
    python3 $sw/VPOT/VPOT.py samplef $out_prefix_other_models $tmp_infile2 $tmp_pedigree
    echo ''

    echo 'Rename the vpot output file to remove the random number in the file name.'
    vpot_output_file_pattern="${out_prefix_other_models}"variant_filtered_output_file_*.txt
    vpot_output_file_renamed="${out_prefix_other_models}"variant_filtered_output_file.txt
    echo 'vpot_output_file=$(ls -t' $vpot_output_file_pattern '| head -n 1)'
    vpot_output_file=$(ls -t $vpot_output_file_pattern | head -n 1)
    echo 'sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7' $vpot_output_file '>' $vpot_output_file_renamed
    sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7 $vpot_output_file > $vpot_output_file_renamed
    echo ''
    rm $vpot_output_file

    echo 'output file:' $vpot_output_file_renamed
    echo ''

else

    # AR autosomal recessive
    # CH compound hets
    # DN denovo

    echo 'Pedigree used for' $inheritance_model':'
    cat $pedigree_file
    echo ''

    echo 'python3' $sw'/VPOT/VPOT.py samplef' $out_prefix_other_models $tmp_infile2 $pedigree_file $proband $inheritance_model
    python3 $sw/VPOT/VPOT.py samplef $out_prefix_other_models $tmp_infile2 $pedigree_file $proband $inheritance_model
    echo ''

    echo 'Rename the vpot output file to remove the random number in the file name.'
    vpot_output_file_pattern="${out_prefix_other_models}${proband}"_"${inheritance_model}"_variant_filtered_output_file_*.txt
    vpot_output_file_renamed="${out_prefix_other_models}${proband}"_"${inheritance_model}"_variant_filtered_output_file.txt
    vpot_output_file=$(ls -t $vpot_output_file_pattern | head -n 1)
    echo 'sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7' $vpot_output_file '>' $vpot_output_file_renamed
    sort -k1,1Vr -k2,2Vr -k3,3V -k4,4V -k6,6 -k7,7 $vpot_output_file > $vpot_output_file_renamed
    echo ''
    rm $vpot_output_file

    echo 'output file:' $vpot_output_file_renamed
    echo ''
fi

echo 'Finished!'

touch "${done_file}"
rm -f "${lock_file}"

