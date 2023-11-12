#!/bin/bash
#PBS -P abcd
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -l mem=4G
#PBS -l ncpus=1
#PBS -l storage=gdata/abcd
set -euo pipefail
: ${chr:? missing variable - chr}
base_d=/my/directory/variant_calling_and_annotation_pipeline/annovar/humandb
d=$base_d/hg19_genome_v211

exit_on_error() {
    if [ $? -ne 0 ] ; then
        echo "TERMINATING.  $1" >&2
        exit 1
    fi
}
cd "$d"
exit_on_error "Unable to change to dir: $d"
#url_pre=https://gnomad-public-us-east-1.s3.amazonaws.com/release/3.1/vcf/genomes
#url_pre=https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/liftover_grch38/vcf/exomes
#url_pre=https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.
url_pre=https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.
url_suf=.vcf.bgz
url_index_suf=.tbi

#f=gnomad.genomes.v3.1.sites.chr${chr}.vcf.bgz
#f=gnomad.exomes.r2.1.1.sites.${chr}.liftover_grch38.vcf.bgz
#wget $url_pre/$f
#exit_on_error "wget ERROR - bgz"
#md5sum "$f" >> /dev/stderr
#wget $url_pre/${f}.tbi
#exit_on_error "wget ERROR - tbi"
f=${url_pre}${chr}$url_suf
wget $f
#md5sum $f >> /dev/stderr
wget ${f}$url_index_suf

echo "SUCCESS" >> /dev/stderr
