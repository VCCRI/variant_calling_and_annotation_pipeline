#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=14:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l storage=gdata/abcd
##PBS -l jobFs=50GB
set -o pipefail
# This program applies some pre-processing to the input
# gnomAD VCF file. This is done to reduce size & retain only 
# relevant information before converting to an annovar database
# (in a future step).
# - values in columns "ID", "QUAL", & "FILTER" are set to ".", as they are not needed
# - only the desired INFO fields are kept
# bcftools is used to do some of the pre-processing & the rest is done via an awk script
#
# Longest job for gnomad 3.1 genomes took ~8hrs
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
#run_ver=$hg38_exomes_211   # completed

scripts_d=/my/directory/variant_calling_and_annotation_pipeline/annovar
out_d=/my/directory/variant_calling_and_annotation_pipeline/annovar/humandb
gnomad_d=/my/directory/variant_calling_and_annotation_pipeline/annovar/humandb
if [ "$run_ver" == "$hg38_exomes_211" ] ; then
    d=$gnomad_d/hg38_exome_v211
    #211 exomes files like:
    #gnomad.exomes.r2.1.1.sites.10.liftover_grch38.vcf.bgz
    in_pre=gnomad.exomes.r2.1.1.sites.
    in_suf=".liftover_grch38.vcf.bgz"
    out_f=$out_d/${in_pre}${chr}.liftover_grch38_processed.vcf
elif [ "$run_ver" == "$hg38_genomes_31" ] ; then
    d=$gnomad_d/hg38_genome_v3.1
    in_pre=gnomad.genomes.v3.1.sites.chr
    in_suf=".vcf.bgz"
    out_f=$out_d/${in_pre}${chr}.processed.vcf
elif [ "$run_ver" == "$hg38_genomes_311" ] ; then
    d=$gnomad_d/hg38_genome_v3.1.1
    in_pre=gnomad.genomes.v3.1.1.sites.chr
    in_suf=".vcf.bgz"
    out_f=$out_d/${in_pre}${chr}.processed.vcf
elif [ "$run_ver" == "$hg19_exomes_211" ] ; then
    d=$gnomad_d/hg19_exome_v211
    #211 exomes files like:
    #gnomad.exomes.r2.1.1.sites.10.vcf.bgz
    in_pre=gnomad.exomes.r2.1.1.sites.
    in_suf=".vcf.bgz"
    out_f=$out_d/${in_pre}${chr}.processed.vcf
fi
cd "$d"
exit_on_error "Unable to change to dir: $d"
input_file=${in_pre}${chr}$in_suf
check_file_exists "$input_file"
# v3.1 genomes
# v2.1.1 exomes
# Remove unwanted information first so that we are working with less data
# Couldn't get bcftools to easily remove unwanted info & set default info tags
# if they didn't exists - so pass through awk script to add default INFO tags
awk_script=$scripts_d/prepare_vcf_file.awk
check_file_exists "$awk_script"
if [ "${run_ver#*_}" == "$exomes_211" ] ; then
    # r211 exomes doesn't store the popmax faf value 
    # - so we will have to calculate it
    faf95_name="faf95_afr,faf95_sas,faf95_amr,faf95_eas,faf95_nfe"
    faf95_default="-999,-999,-999,-999,-999"
else
    # for r31 genomes:
    faf95_name=faf95_popmax
    faf_95_default="-999"
fi
# write field names & defaults to a variable & pass to awk script - as name can
# change between gnomad versions.  The variable contents will be in order of 
# desired output
info_fields_str=AC,AF,AF_popmax,popmax,nhomalt,$faf95_name
# store variable in required format for bcftools command:
#INFO/AC,INFO/AF,INFO/AF_popmax,INFO/popmax,INFO/nhomalt,INFO/faf95
bcf_info_fields_str=$( awk -F"," '{for(i=1; i<NF;++i) printf "INFO/%s,", $i ; printf "INFO/%s", $NF}' <<< "$info_fields_str" )
info_defaults_str="-999,-999,-999,.,-999,$faf95_default"
# use bcftools annotate because it tidies up the header, and is
# simple to use.  
bcftools annotate "${input_file}" -x ID,QUAL,FILTER,^$bcf_info_fields_str \
    | awk -f "$awk_script" -v info_fields_str="$info_fields_str" \
    -v info_defaults_str="$info_defaults_str" > "$out_f"
exit_on_error "ERROR processing file"
echo "SUCCESS" >>/dev/stderr
