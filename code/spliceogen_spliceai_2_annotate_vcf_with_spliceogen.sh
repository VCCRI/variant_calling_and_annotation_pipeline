#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=96G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N annotatespliceogen
#PBS -lstorage=gdata/abcd

#input parameters:
#sample=$1
#infile1=$2
#infile2=$3
#outdir=$4
#outfile=$5
#sw_and_refs=$6

set -euo pipefail

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"
tmpdir_for_bcftools="${tmpdir}"/tmp_for_bcftools
mkdir -p "${tmpdir_for_bcftools}"
tmpdir_for_annovar="${tmpdir}"/humandb
mkdir -p "${tmpdir_for_annovar}"

queue_file="${outfile}.queued"
lock_file="${outfile}.lock"
done_file="${outfile}.done"
term_file="${outfile}.term"
log_file="${outfile}.log"

term_handler()
{
#    rm -f "${lock_file}"
    touch "${term_file}"
    exit 1
}
trap 'term_handler' TERM

touch "${lock_file}"
rm -f "${queue_file}"

module load python3
module load bcftools

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmp_spliceogen_annovar_table="${tmpdir}"/tmp_spliceogen_annovar_table.vcf
tmp_spliceogen_annovar_table_sorted="${tmpdir}"/tmp_spliceogen_annovar_table_sorted.vcf
tmp_spliceogen_annovar_table_sorted_convert2annovar_avinput="${tmpdir}"/tmp_spliceogen_annovar_table_sorted_convert2annovar.avinput
tmp_spliceogen_annovar_table_sorted_txt="${tmpdir_for_annovar}"/hg38_spliceogen.txt
infile1_basename=$(basename $infile1)
tmp_infile1_gz="${tmpdir}"/"${infile1_basename}"
tmp_infile1="${tmp_infile1_gz%.gz}"

echo ''
echo '##### Convert' $infile2 'to' $tmp_spliceogen_annovar_table_sorted_txt 'so that it can be input to annovar to turn it into an annovar reference table'
echo ''

echo 'Build header for' $tmp_spliceogen_annovar_table
echo -e '##fileformat=VCFv4.2' > $tmp_spliceogen_annovar_table
echo -e '##fileDate=2021-01-01' >> $tmp_spliceogen_annovar_table
echo -e '##source=Spliceogen' >> $tmp_spliceogen_annovar_table
echo -e '##reference=GRCh38' >> $tmp_spliceogen_annovar_table
echo -e '##INFO=<ID=spliceogen_gene,Number=.,Type=String,Description="Spliceogen gene">' >> $tmp_spliceogen_annovar_table
echo -e '##INFO=<ID=spliceogen_withinSite,Number=.,Type=String,Description="Spliceogen withinSite">' >> $tmp_spliceogen_annovar_table
echo -e '##INFO=<ID=spliceogen_donGainP,Number=.,Type=String,Description="Spliceogen donGainP">' >> $tmp_spliceogen_annovar_table
echo -e '##INFO=<ID=spliceogen_accGainP,Number=.,Type=String,Description="Spliceogen accGainP">' >> $tmp_spliceogen_annovar_table
echo -e '##INFO=<ID=spliceogen_donLossP,Number=.,Type=String,Description="Spliceogen donLossP">' >> $tmp_spliceogen_annovar_table
echo -e '##INFO=<ID=spliceogen_accLossP,Number=.,Type=String,Description="Spliceogen accLossP">' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr1,length=248956422>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr2,length=242193529>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr3,length=198295559>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr4,length=190214555>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr5,length=181538259>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr6,length=170805979>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr7,length=159345973>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr8,length=145138636>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr9,length=138394717>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr10,length=133797422>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr11,length=135086622>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr12,length=133275309>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr13,length=114364328>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr14,length=107043718>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr15,length=101991189>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr16,length=90338345>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr17,length=83257441>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr18,length=80373285>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr19,length=58617616>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr20,length=64444167>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr21,length=46709983>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chr22,length=50818468>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chrX,length=156040895>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chrY,length=57227415>' >> $tmp_spliceogen_annovar_table
echo -e '##contig=<ID=chrM,length=16569>' >> $tmp_spliceogen_annovar_table
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> $tmp_spliceogen_annovar_table
echo ''

count=`grep -v '^#' $infile2 | grep -v '^CHROM' | grep -v '^chrom' | wc -l | cut -d' ' -f1 || true`
echo $count 'data records in' $infile2
echo ''

col_gene=0
col_withinSite=0
col_donGainP=0
col_accGainP=0
col_donLossP=0
col_accLossP=0
echo 'col_gene_result=head -n 1' $infile2 '| sed -e s/\t/\n/g | grep -n GENE'
col_gene_result=`head -n 1 $infile2 | sed -e 's/\t/\n/g' | grep -n GENE`
echo ''
echo 'col_withinSite_result=head -n 1' $infile2 '| sed -e s/\t/\n/g | grep -n withinSite'
col_withinSite_result=`head -n 1 $infile2 | sed -e 's/\t/\n/g' | grep -n withinSite`
echo ''
echo 'col_donGainP_result=head -n '1 $infile2 '| sed -e s/\t/\n/g | grep -n donGainP'
col_donGainP_result=`head -n 1 $infile2 | sed -e 's/\t/\n/g' | grep -n donGainP`
echo ''
echo 'col_accGainP_result=head -n '1 $infile2 '| sed -e s/\t/\n/g | grep -n accGainP'
col_accGainP_result=`head -n 1 $infile2 | sed -e 's/\t/\n/g' | grep -n accGainP`
echo ''
echo 'col_donLossP_result=head -n '1 $infile2 '| sed -e s/\t/\n/g | grep -n donLossP'
col_donLossP_result=`head -n 1 $infile2 | sed -e 's/\t/\n/g' | grep -n donLossP`
echo ''
echo 'col_accLossP_result=head -n '1 $infile2 '| sed -e s/\t/\n/g | grep -n accLossP'
col_accLossP_result=`head -n 1 $infile2 | sed -e 's/\t/\n/g' | grep -n accLossP`
echo ''
IFS=':' read -r -a array <<< "$col_gene_result"
col_gene="${array[0]}"
IFS=':' read -r -a array <<< "$col_withinSite_result"
col_withinSite="${array[0]}"
IFS=':' read -r -a array <<< "$col_donGainP_result"
col_donGainP="${array[0]}"
IFS=':' read -r -a array <<< "$col_accGainP_result"
col_accGainP="${array[0]}"
IFS=':' read -r -a array <<< "$col_donLossP_result"
col_donLossP="${array[0]}"
IFS=':' read -r -a array <<< "$col_accLossP_result"
col_accLossP="${array[0]}"
echo 'col_gene =' $col_gene
echo 'col_withinSite =' $col_withinSite
echo 'col_donGainP =' $col_donGainP
echo 'col_accGainP =' $col_accGainP
echo 'col_donLossP =' $col_donLossP
echo 'col_accLossP =' $col_accLossP
echo ''

#echo 'grep -v' '^#' $infile2 '| sed' 's/;/|/g' '| awk -v col_gene='$col_gene '-v col_withinSite='$col_withinSite '-v col_donGainP='$col_donGainP '-v col_accGainP='$col_accGainP '-v col_donLossP='$col_donLossP '-v col_accLossP='$col_accLossP 'BEGIN {FS="\t";OFS="\t"} {if ( (($col_withinSite!=".") && ($col_withinSite!="")) && (($col_donGainP>0.7) || ($col_accGainP>0.7) || ($donLossP>0.8) || ($accLossP>0.8)) ) {print $1, $2, ".", $4, $5, ".", ".", "spliceogen_gene=" $col_gene ";spliceogen_withinSite=" $withinSite ";spliceogen_donGainP=" $donGainP ";spliceogen_accGainP=" $col_accGainP ";spliceogen_donLossP=" $col_donLossP ";spliceogen_accLossP=" $col_accLossP}}' '>>' $tmp_spliceogen_annovar_table
#grep -v '^#' $infile2 | sed 's/;/|/g' | awk -v col_gene=$col_gene -v col_withinSite=$col_withinSite -v col_donGainP=$col_donGainP -v col_accGainP=$col_accGainP -v col_donLossP=$col_donLossP -v col_accLossP=$col_accLossP 'BEGIN {FS="\t";OFS="\t"} {if ( (($col_withinSite!=".") && ($col_withinSite!="")) && (($col_donGainP>0.7) || ($col_accGainP>0.7) || ($col_donLossP>0.8) || ($col_accLossP>0.8)) ) {print $1, $2, ".", $4, $5, ".", ".", "spliceogen_gene=" $col_gene ";spliceogen_withinSite=" $col_withinSite ";spliceogen_donGainP=" $col_donGainP ";spliceogen_accGainP=" $col_accGainP ";spliceogen_donLossP=" $col_donLossP ";spliceogen_accLossP=" $col_accLossP}}' >> $tmp_spliceogen_annovar_table
#echo ''

echo 'grep -v' '^#' $infile2 '| sed' 's/;/|/g' '| awk -v col_gene='$col_gene '-v col_withinSite='$col_withinSite '-v col_donGainP='$col_donGainP '-v col_accGainP='$col_accGainP '-v col_donLossP='$col_donLossP '-v col_accLossP='$col_accLossP 'BEGIN {FS="\t";OFS="\t"} {print $1, $2, ".", $4, $5, ".", ".", "spliceogen_gene=" $col_gene ";spliceogen_withinSite=" $col_withinSite ";spliceogen_donGainP=" $col_donGainP ";spliceogen_accGainP=" $col_accGainP ";spliceogen_donLossP=" $col_donLossP ";spliceogen_accLossP=" $col_accLossP}' '>>' $tmp_spliceogen_annovar_table
grep -v '^#' $infile2 | sed 's/;/|/g' | awk -v col_gene=$col_gene -v col_withinSite=$col_withinSite -v col_donGainP=$col_donGainP -v col_accGainP=$col_accGainP -v col_donLossP=$col_donLossP -v col_accLossP=$col_accLossP 'BEGIN {FS="\t";OFS="\t"} {print $1, $2, ".", $4, $5, ".", ".", "spliceogen_gene=" $col_gene ";spliceogen_withinSite=" $col_withinSite ";spliceogen_donGainP=" $col_donGainP ";spliceogen_accGainP=" $col_accGainP ";spliceogen_donLossP=" $col_donLossP ";spliceogen_accLossP=" $col_accLossP}' >> $tmp_spliceogen_annovar_table
echo ''

echo 'bcftools sort -o' $tmp_spliceogen_annovar_table_sorted $tmp_spliceogen_annovar_table '-T' $tmpdir_for_bcftools
bcftools sort -o $tmp_spliceogen_annovar_table_sorted $tmp_spliceogen_annovar_table -T $tmpdir_for_bcftools
echo ''

count=`grep -v '^#' $tmp_spliceogen_annovar_table_sorted | grep -v '^CHROM' | grep -v '^chrom' | wc -l | cut -d' ' -f1 || true`
echo $count 'data records in' $tmp_spliceogen_annovar_table_sorted
echo ''

echo 'perl' "${sw}"'/annovar/convert2annovar.pl -format vcf4' $tmp_spliceogen_annovar_table_sorted '>' $tmp_spliceogen_annovar_table_sorted_convert2annovar_avinput
perl "${sw}"/annovar/convert2annovar.pl -format vcf4 $tmp_spliceogen_annovar_table_sorted > $tmp_spliceogen_annovar_table_sorted_convert2annovar_avinput
echo ''

echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_and_its_annovar_avinput_files_to_tab_delimited_for_annovar.py -i' $tmp_spliceogen_annovar_table_sorted '-a' $tmp_spliceogen_annovar_table_sorted_convert2annovar_avinput '-o' $tmp_spliceogen_annovar_table_sorted_txt
python3 "${sw}"/victorchang_scripts/convert_vcf_and_its_annovar_avinput_files_to_tab_delimited_for_annovar.py -i $tmp_spliceogen_annovar_table_sorted -a $tmp_spliceogen_annovar_table_sorted_convert2annovar_avinput -o $tmp_spliceogen_annovar_table_sorted_txt
echo ''

echo 'perl' "${sw}"'/annovar/index_annovar.pl' $tmp_spliceogen_annovar_table_sorted_txt '-outfile' "${tmp_spliceogen_annovar_table_sorted_txt}"'.idx --skipsort'
perl "${sw}"/annovar/index_annovar.pl $tmp_spliceogen_annovar_table_sorted_txt -outfile "${tmp_spliceogen_annovar_table_sorted_txt}".idx --skipsort
echo ''

echo 'mv' "${tmp_spliceogen_annovar_table_sorted_txt}".idx.idx "${tmp_spliceogen_annovar_table_sorted_txt}"'.idx'
mv "${tmp_spliceogen_annovar_table_sorted_txt}".idx.idx "${tmp_spliceogen_annovar_table_sorted_txt}".idx
echo ''

echo ''
echo '##### Annotate vcf file' $infile1 'with its spliceogen values in table' $tmp_spliceogen_annovar_table_sorted_txt 'to produce' $outfile 'so that the vcf contains spliceogen annotations in addition to other annotations'
echo ''

if [[ $infile1 == *.vcf.gz ]]; then
  #echo 'cp' $infile1 $tmp_infile1_gz
  #cp $infile1 $tmp_infile1_gz
  #echo ''
  #echo 'gunzip -f' $tmp_infile1_gz
  #gunzip -f $tmp_infile1_gz
  #echo ''
  tmp_infile1=$infile1
else
  tmp_infile1=$infile1
fi

echo 'perl' $sw'/annovar/table_annovar.pl' $tmp_infile1 "${tmpdir_for_annovar}"'/ -vcfinput -buildver' $genome_version '\'
echo '  -out' $outfile '-remove \'
echo '  -protocol' spliceogen '\'
echo '  -operation f -nastring . \'
echo '  -arg -time'
perl $sw/annovar/table_annovar.pl $tmp_infile1 "${tmpdir_for_annovar}"/ -vcfinput -buildver "$genome_version" \
  -out $outfile -remove \
  -protocol spliceogen \
  -operation f -nastring . \
  -arg -time
echo ''

echo 'python3' "${sw}"'/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i' "${outfile}"'."${genome_version}"_multianno.vcf -o' "${outfile}"'.tsv -end_id YES'
python3 "${sw}"/victorchang_scripts/convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i "${outfile}"."${genome_version}"_multianno.vcf -o "${outfile}".tsv -end_id YES
echo ''

echo ''
echo 'outfile1:' "${outfile}"."${genome_version}"_multianno.vcf
echo 'outfile2:' "${outfile}".tsv
echo ''
echo 'Finished!'
echo ''



