#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=6:00:00
#PBS -l mem=11G
#PBS -l ncpus=1
#PBS -l jobfs=10G
#PBS -lstorage=scratch/abcd+gdata/abcd
#noPBS -m bea

# memory is usually 8G
# walltime is usually 6hrs

exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "Terminating.  $1" >&2
        exit 1
    fi
}

# GT 1/. causes plink --missing to crash, so change to 0/1

sw_and_refs=/my/directory/variant_calling_and_annotation_pipeline/code/where_are_softwares_and_references.sh
. "${sw_and_refs}"

sample=$batch
infile1=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".AllChrom.GT01.vcf.gz # This is original input vcf, after converting GT=1/. to GT=0/1, before vcfQC_1_copy converted it to VCFv4.0
infile2=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".AllChrom.GT01_VCFv4p0.vcf.gz # This is input vcf converted it to VCFv4.0 and also GT=1/. converted to GT=0/1
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc
outfile=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/"${batch}".plink_output

# directories & filenames
vcf_d=$outdir
in_vcf=$infile1
in_vcf_4_0=$infile2

#queue_file="${outfile}.queued"
#lock_file="${outfile}.lock"
#done_file="${outfile}.done"
#term_file="${outfile}.term"
#log_file="${outfile}.log"

#touch "${lock_file}"
#rm -f "${queue_file}"

module load java/jdk-8.40
module load samtools

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

##### script directory
fam_ped_script=/my/directory/variant_calling_and_annotation_pipeline/victorchang_scripts/vcfQC_modify_fam_ped.awk
plink_ped_script=/my/directory/variant_calling_and_annotation_pipeline/victorchang_scripts/vcfQC_modify_plink_ped.awk

##### software
pseq=$sw/plinkseq-0.10/pseq
plink=$sw/plink-1.07-x86_64/plink
king=$sw/KING/2.1.4/king
vcftools=$sw/vcftools-0.14/bin/vcftools

##### ped file:
#     Family ID
#     Individual ID
#     Paternal ID
#     Maternal ID
#     Sex (1=male; 2=female; other=unknown)
#     Phenotype
echo 'create a ped file for this singleton input'
ped_file="${outdir}"/"${sample}".ped
ped_file_master=/my/directory/variant_calling_and_annotation_pipeline/code/"${batch}".pedigree_file.ped
cp $ped_file_master $ped_file

##### set up output files from vcftools --plink
echo 'set up output files from vcftools --plink'
out_d="${outdir}"/"${sample}".plink
mkdir -p "${out_d}"
cd "$out_d"
exit_on_error "Unable to change to output directory: $out_d"
ped=$ped_file
ls "${ped}" &>/dev/null
exit_on_error "Unable to locate ped file: ${ped}"
chrall="${out_d}"/"${sample}".chrall
fam=${chrall}.fam
original_fam=${chrall}.fam.original
plink_ped=${chrall}.ped
original_plink_ped=${plink_ped}.original
original_fam=${fam}.original
bed=${chrall}.bed
missing=$out_d/"${sample}".chrall_missing
chksex=$out_d/"${sample}".chrall_chksex
ls "${in_vcf_4_0}" &>/dev/null
exit_on_error "Unable to locate in_vcf_4_0: ${in_vcf_4_0}"
ls "$in_vcf" &>/dev/null
exit_on_error "Unable to locate in_vcf: $in_vcf"
echo ''

# plink istats
echo 'plink istats'
echo $pseq "${in_vcf_4_0}" 'i-stats >' "${out_d}"'/'"${sample}"'_vcf_istats.txt'
$pseq "${in_vcf_4_0}" i-stats > "${out_d}"/"${sample}"_vcf_istats.txt
echo ''

# vcftools plink
echo 'vcftools plink'
echo $vcftools '--gzvcf' $in_vcf '--plink --temp' $tmpdir '--out' $chrall
$vcftools --gzvcf $in_vcf --plink --temp $tmpdir --out $chrall
exit_on_error "vcftools --plink unsuccessful."
echo ''

# vcftools --plink creates a .fam file - we need to fill in missing ped details
# NOTE that the .fam file generally gets overwritten from the information found
# in the plink ped file ; just doing this modification as a precaution
# check that the fam file exists
# This file is probably not created for samples whose pedigree say that they are not part of a family.
if [[ -f "$fam" ]]; then

  #echo 'vcftools --plink creates a .fam file - we need to fill in missing ped details'
  #ls "$fam" &>/dev/null
  #exit_on_error "Unable to locate fam: $fam"
  #echo ''

  # make a copy of the original fam file
  echo 'make a copy of the original fam file'
  cp "$fam" "$original_fam"
  exit_on_error "Error copying fam file to .original"
  echo 'awk -v fam='$original_fam '-v ped='$ped '-f' $fam_ped_script '>' $fam
  awk -v fam=$original_fam -v ped=$ped -f $fam_ped_script > $fam
  exit_on_error "modfiy fam ped unsuccessful"
  echo ''
fi

# vcftools --plink creates a .ped file - we need to fill in missing ped details
echo 'vcftools --plink creates a .ped file - we need to fill in missing ped details'
# check that the plink ped file exists
ls "$plink_ped" &>/dev/null
exit_on_error "Unable to locate plink_ped: $plink_ped"
echo ''

# make a copy of the original plink_ped file
echo 'make a copy of the original plink_ped file'
cp "$plink_ped" "$original_plink_ped"
exit_on_error "Error copying plink_ped file to .original"
echo 'awk -v plink_ped='$original_plink_ped '-v our_ped='$ped '-f' $plink_ped_script '>' $plink_ped
awk -v plink_ped=$original_plink_ped -v our_ped=$ped -f $plink_ped_script > $plink_ped
exit_on_error "modfiy plink ped unsuccessful"
echo ''

# plink missing # this step crashes if input GT = 1/.
echo $plink '--file' $chrall '--noweb --allow-no-sex --missing --out' $missing
$plink --file $chrall --noweb --allow-no-sex --missing --out $missing
echo ''

# make bed ; NOTE that this command overwrites the .fam file
echo $plink '--file' $chrall '--noweb --allow-no-sex --make-bed --out' $chrall
$plink --file $chrall --noweb --allow-no-sex --make-bed --out $chrall
echo ''

# check sex
echo $plink '--bfile' $chrall '--noweb --allow-no-sex --check-sex --out' $chksex
$plink --bfile $chrall --noweb --allow-no-sex --check-sex --out $chksex
echo ''

# kinship check using KING
echo $king '-b' $bed '--kinship'
$king -b $bed --kinship
echo ''
echo $king '-b' $bed '--ibdseg'
$king -b $bed --ibdseg
echo ''
echo $king '-b' $bed '--ibs'
$king -b $bed --ibs
echo ''

echo ''
echo 'Finished!'
echo ''

#touch "${done_file}"
#rm -f "${lock_file}"



