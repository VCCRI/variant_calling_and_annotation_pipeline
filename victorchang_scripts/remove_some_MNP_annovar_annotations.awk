# Author: Emma Rath
# Date: 2021-01-21
#
# awk -f remove_some_MNP_annovar_annotations.awk input.vcf > output.vcf
#
# For multinucleotide polymorphisms/variants (MNPs/MNVs) not in normalised form,
# the derived annovar fields of consequences are incorrect.
# Thus, remove those annotations.
# MNPs have length of REF field > 1 and length of ALT field > 1.
#
# Here are some examples of incorrect annotations:
# hg38 chr20:62324148-62324150:CAC>CAT is annotated with c.5698delG:p.V1900X stopgain
# hg38 chr7:92469384-92469386:AGT>CAT is annotated with c.997_998TG
# hg38 chr7:92469602-92469604:TTG>TCA is annotated with c.778_779TG
# hg38 chr17:46033173-46033175:GCG>GCA is annotated with c.2550delC:p.D850fs frameshift_deletion
#
# The fields to remove are:
# Func.refGene				exonic
# GeneDetail.refGene			dist\x3d898792,dist\x3d375609
# ExonicFunc.refGene			stopgain
# AAChange.refGene			LAMA5:NM_005560:exon43:c.5698delG:p.V1900X
# Func.wgEncodeGencodeBasicV33		exonic
# GeneDetail.wgEncodeGencodeBasicV33	ENST00000551531.1:exon3:c.177-2A>T,ENST00000551940.1:exon4:c.344-2A>T
# ExonicFunc.wgEncodeGencodeBasicV33	stopgain
# AAChange.wgEncodeGencodeBasicV33	LAMA5:ENST00000252999.7:exon43:c.5698delG:p.V1900X
# VARIANT_TYPE				exonic_stopgain_
# For the GencodeBasic field, remove regardless of version
# The refGene fields are created when refGene annovar annoation is used.
# The wgEncodeGencodeBasicV33 fields are created when Gencode annovar annotation is used.
# The VARIANT_TYPE field is created by table_annovar_ed.pl, not by table_annovar.pl

function fix_info(old_info) {
  new_info = ""
  num_in_array = split(old_info, array, ";")
  for (j=1; j<=num_in_array; ++j) {
    this_field = array[j]
    rmv_field = 0
    if (substr(this_field,1,12)=="Func.refGene") { rmv_field = 1 }
    if (substr(this_field,1,18)=="GeneDetail.refGene") { rmv_field = 1 }
    if (substr(this_field,1,18)=="ExonicFunc.refGene") { rmv_field = 1 }
    if (substr(this_field,1,16)=="AAChange.refGene") { rmv_field = 1 }
    if (substr(this_field,1,25)=="Func.wgEncodeGencodeBasic") { rmv_field = 1 }
    if (substr(this_field,1,31)=="GeneDetail.wgEncodeGencodeBasic") { rmv_field = 1 }
    if (substr(this_field,1,31)=="ExonicFunc.wgEncodeGencodeBasic") { rmv_field = 1 }
    if (substr(this_field,1,29)=="AAChange.wgEncodeGencodeBasic") { rmv_field = 1 }
    if (substr(this_field,1,12)=="VARIANT_TYPE") { rmv_field = 1 }
    if (rmv_field == 1) {
      split(this_field, array2, "=")
      this_field = array2[1]"=."
    }
    if (new_info == "") {
      new_info = this_field
    } else {
      new_info = new_info";"this_field
    }
  }
  return new_info
}

BEGIN {
    FS = "\t"
    OFS = "\t"
}
{
  char1 = substr($1,1,1)
  if (char1 == "#") {
    print $0
  } else {
    if ((length($4)>1) && (length($5)>1)) {
      new_info = fix_info($8)
      $8 = new_info
      print $0
    } else {
      print $0
    }
  }
}


