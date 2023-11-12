#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=04:00:00
#PBS -l mem=4GB
#PBS -l ncpus=1
#PBS -l storage=gdata/abcd
#PBS -l jobFs=64GB
set -uo pipefail
# This file concatinates all the chromosome annovar db's
# , checks they have the same headers, and creates annovar index
#
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
last_5_cols_faf95_pop=0
scripts_d=/my/directory/variant_calling_and_annotation_pipeline/annovar
d=/my/directory/variant_calling_and_annotation_pipeline/annovar/humandb
cd "$d"
exit_on_error "Unable to change to dir: $d"
sw=/g/data/jb96/software
index_pl=$sw/annovar/index_annovar.pl
check_file_exists "$index_pl"
# prefix for annovar db fields #TODO - change as appropriate
################
# gnomad 3.1 genome
################
#db_pre=gnomad31_genome_VC_
################
################
# gnomad 2.1.1 exome
################
#db_pre=gnomad211_exome_VC_
################
# input file like:
#TODO - change as appropriate
################
# gnomad 3.1 genome
################
# gnomad.genomes.v3.1.sites.chrY.processed.vcf
#in_file=gnomad.genomes.v3.1.sites.chr${chr}.processed.vcf
in_pre="gnomad.genomes.v3.1.sites.chr"
in_suf=".processed.vcf_annovar.txt"
out_f=hg38_gnomad31_VC_genome.txt
################
################
# gnomad 2.1.1 exome
################
# gnomad.exomes.r2.1.1.sites.7.liftover_grch38_processed.vcf
#in_file=gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38_processed.vcf
#in_pre="gnomad.exomes.r2.1.1.sites."
#for hg38:
#in_suf=".liftover_grch38_processed.vcf_annovar.txt"
#out_f=hg38_gnomad211_VC_exome.txt
#for hg19:
#in_suf=".processed.vcf_annovar.txt"
#out_f=hg19_gnomad211_VC_exome.txt
# gnomad doesn't supply faf95_popmax for v211 - we have to caluculate
# it from the individual populations
#last_5_cols_faf95_pop=1
################

# create some temporary filenames
tmp_out_f=tmp_$out_f
tmp_out_post_processed_f=tmp_post_processed_$out_f

# first check that all the headers are the same
chr=1
in_file=${in_pre}${chr}${in_suf}
check_file_exists "$in_file"
head -1 "$in_file" > tmp_hdr_$chr
{
    cat "$in_file"
    for chr in {2..22} X Y ; do
        in_file=${in_pre}${chr}${in_suf}
        check_file_exists "$in_file"
        head -1 "$in_file" > tmp_hdr_$chr
        diff -q tmp_hdr_1 tmp_hdr_$chr &>/dev/null
        if [ $? -ne 0 ] ; then
            rm tmp_hdr_*
            echo "ERROR.  mismatched header for chr $chr" >>/dev/stderr
            exit 1
        fi
        tail -n+2 "$in_file"
        exit_on_error "Problem outputing for $in_file"
    done
} > "$tmp_out_f"
if [ $last_5_cols_faf95_pop -eq 1 ] ; then
    #TODO
    # deal with header & data - calculate faf95_popmax
    # from individual population data
    awk 'BEGIN {
        FS = "\t"
        OFS = "\t"
    }
    # header
    NR==1 {
        for (i=1; i<=NF-5; ++i)
            printf "%s\t", $i
        # col names like: gnomad211_exome_faf95_nfe
        printf "%s\n", "gnomad211_exome_faf95_popmax"
        next
    }
    {
        # calc popmax
        popmax = -999
        for (i=NF-5+1; i<=NF; ++i) {
            if ($i > popmax)
                popmax = $i
        }
        for (i=1; i<=NF-5; ++i)
            printf "%s\t", $i
        printf "%s\n", popmax
    }' "$tmp_out_f" > "$tmp_out_post_processed_f"
else
    tmp_out_post_processed_f=$tmp_out_f
fi

# tidy up
rm tmp_hdr_*
# index file
echo "About to index.." >> /dev/stderr
perl "$index_pl" "$tmp_out_post_processed_f" -outfile "$out_f"
exit_on_error "Problem indexing $out_f"
rm "$tmp_out_f"
rm "$tmp_out_post_processed_f"
echo "SUCCESS" >>/dev/stderr
