#!/bin/bash
#PBS -P abcd
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=32G
#PBS -l ncpus=1
#PBS -l jobfs=1G
#PBS -N annovar
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

# set environment variables for softwares and references
# sw_and_refs=where_are_softwares_and_references.sh
. "${sw_and_refs}"

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

echo ''
echo "Start running annovar for" $infile
echo ''

# Most of these annotations come from annovar
# described at https://doc-openbio.readthedocs.io/projects/annovar/en/latest/user-guide/download/
# Others, such as BayesDel and latest clinvar, have been created manually inhouse for annovar.

if [[ $genome_version == "hg38" ]]; then

  # annotation cmd with all - for production usage
  echo $sw'/annovar/table_annovar_ed.pl' "$infile" ${humandb}'/ -vcfinput -buildver' $genome_version '\'
  echo '  -out' "$outfile" '-remove \'
  echo '  -protocol refGene,'$wgEncodeGencodeBasic','$EHRF_regions',cytoBand,genomicSuperDups,FANTOM5_CAGEr_humanHeart_transcriptStartSites,hotspot_exons,gwava,cadd_v1p6,gnomad211_VC_exome,gnomad31_VC_genome,gnomadMNV,gnomadMNVdist,'$exac03_ed',exac03_HOMO,WGSb1_5_full_freq_final_v0,esp6500siv2_all,UK10K,mcap,revel,LINSIGHT_ed,clinvar_20221008,BayesDel_addAF_noAF,pext,avsnp150,dbnsfp42a_minus_some_cols_2021aug,dbscsnv11,esp6500siv2_aa,esp6500siv2_ea,gene4denovo201907,gme,hrcr1,intervar_20180118,kaviar_20150923,regsnpintron,mutpred_indels,spliceai_raw_noZeros,spliceaiIndels_raw_noZeros,spliceai_masked_noZeros,spliceaiIndels_masked_noZeros,WGSb1_8_full_freq,mmsplice_ge2,SynVarDisruptMrna_SURF_ge_5_or_synvar,EVE \'
  echo '  -operation g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . \'
  echo '  -arg -exonicsplicing,-exonicsplicing,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time'

  $sw/annovar/table_annovar_ed.pl "$infile" ${humandb}/ -vcfinput -buildver "$genome_version" \
    -out "$outfile" -remove \
    -protocol refGene,"$wgEncodeGencodeBasic","$EHRF_regions",cytoBand,genomicSuperDups,FANTOM5_CAGEr_humanHeart_transcriptStartSites,hotspot_exons,gwava,cadd_v1p6,gnomad211_VC_exome,gnomad31_VC_genome,gnomadMNV,gnomadMNVdist,"$exac03_ed",exac03_HOMO,WGSb1_5_full_freq_final_v0,esp6500siv2_all,UK10K,mcap,revel,LINSIGHT_ed,clinvar_20221008,BayesDel_addAF_noAF,pext,avsnp150,dbnsfp42a_minus_some_cols_2021aug,dbscsnv11,esp6500siv2_aa,esp6500siv2_ea,gene4denovo201907,gme,hrcr1,intervar_20180118,kaviar_20150923,regsnpintron,mutpred_indels,spliceai_raw_noZeros,spliceaiIndels_raw_noZeros,spliceai_masked_noZeros,spliceaiIndels_masked_noZeros,WGSb1_8_full_freq,mmsplice_ge2,SynVarDisruptMrna_SURF_ge_5_or_synvar,EVE \
    -operation g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . \
    -arg -exonicsplicing,-exonicsplicing,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time

fi

if [[ $genome_version == "hg19" ]]; then

  # annotation cmd with all - for production usage
  echo $sw'/annovar/table_annovar_ed.pl' "$infile" ${humandb}'/ -vcfinput -buildver' $genome_version '\'
  echo '  -out' "$outfile" '-remove \'
  echo '  -protocol refGene,'$wgEncodeGencodeBasic','$EHRF_regions',cytoBand,genomicSuperDups,FANTOM5_CAGEr_humanHeart_transcriptStartSites,hotspot_exons,gwava,gnomad211_VC_exome,gnomad31_VC_genome,gnomadMNV,gnomadMNVdist,'$exac03_ed',esp6500siv2_all,revel,clinvar_20221008,BayesDel_all,pext,avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,regsnpintron,spliceai_raw_noZeros,spliceaiIndels_raw_noZeros,spliceai_masked_noZeros,spliceaiIndels_masked_noZeros,mmsplice_ge2,SynVarDisruptMrna_SURF_ge_5_or_synvar \'
  echo '  -operation g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . \'
  echo '  -arg -exonicsplicing,-exonicsplicing,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time'

  $sw/annovar/table_annovar_ed.pl "$infile" ${humandb}/ -vcfinput -buildver "$genome_version" \
    -out "$outfile" -remove \
    -protocol refGene,"$wgEncodeGencodeBasic","$EHRF_regions",cytoBand,genomicSuperDups,FANTOM5_CAGEr_humanHeart_transcriptStartSites,hotspot_exons,gwava,gnomad211_VC_exome,gnomad31_VC_genome,gnomadMNV,gnomadMNVdist,"$exac03_ed",esp6500siv2_all,revel,clinvar_20221008,BayesDel_all,pext,avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,regsnpintron,spliceai_raw_noZeros,spliceaiIndels_raw_noZeros,spliceai_masked_noZeros,spliceaiIndels_masked_noZeros,mmsplice_ge2,SynVarDisruptMrna_SURF_ge_5_or_synvar \
    -operation g,g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . \
    -arg -exonicsplicing,-exonicsplicing,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time,-time

fi

# actual outfile from annovar is "${outfile}"."${genome_version}"_multianno.vcf
outfile1="${outfile}"."${genome_version}"_multianno.vcf
outfile2="${outfile}"."${genome_version}"_multianno.rmvSomeAnnovar.vcf

echo 'awk -f' $sw'/victorchang_scripts/remove_some_MNP_annovar_annotations.awk' $outfile1 '>' $outfile2
awk -f $sw/victorchang_scripts/remove_some_MNP_annovar_annotations.awk $outfile1 > $outfile2
echo ''

touch "${done_file}"
rm -f "${lock_file}"

echo ''
echo 'Finished!'
echo ''
