# Author: Emma Rath
# Date: 2020-12-09
#
# awk -f reformat_vcf_for_vpot_and_set_NR_NV_to_max_values.awk temp.vcf > temp_out.vcf
#
# The sample FORMAT field GT can contain slash or bar (eg. 0/1 or 0|1).
# Vpot needs it to be slash (eg. 0/1) and will crash if it is bar (eg. 0|1).
# Thus, change the bar to slash.
#
# Here is an explanation of the fix for platypus NR and NV
#
# Platypus does not fill in the definition of FORMAT fields.
# when multiple alleles are split into multiple lines of single alleles, the record counts/depths don't get transferred,
# and eventually these variants get a GT of dot instead of 1 for 0/1 in VPOT, and so the variants are lost.
#
# Platypus headers have the following:
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=GQ,Number=.,Type=Integer,Description="Genotype quality as phred score">
##FORMAT=<ID=GOF,Number=.,Type=Float,Description="Goodness of fit value">
##FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">
##FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
#
# This means that NR and NV fields do not get their values split by vt decompose,
# so that multiple allelic values remain even when the alt has only one value (vt decompose converts multi-allelic to multiple single-allele rows)
# This would cause problems downstream in vpot when vpot needs to use NR and NV to calculate depth.
#
# Instead they should be the following:
##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotypes">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality as phred score">
##FORMAT=<ID=GOF,Number=1,Type=Float,Description="Goodness of fit value">
##FORMAT=<ID=NR,Number=A,Type=Integer,Description="Number of reads covering variant location in this sample">
##FORMAT=<ID=GL,Number=R,Type=Float,Description="Genotype log10-likelihoods for AA,AB and BB genotypes, where A = ref and B = variant. Only applicable for bi-allelic sites">
##FORMAT=<ID=NV,Number=A,Type=Integer,Description="Number of reads containing variant in this sample">
#
# To fix this, this program looks to see if there is an NR or NV sample field.
# If there is, and if is has multiple values separated by a comma,
# then this program puts only one value in the field - the highest value.

function find_NR_NV_col(data, col_name, delim) {
  col_num = 0
  num_cols = split(data, arr, delim)
  for (z=1; z<=num_cols; z++) {
    this_col_name=arr[z]
    if (this_col_name==col_name) {
      col_num = z
    }
  }
  return col_num
}

function fix_GT_col(data, col_num, delim) {
  num_cols = split(data, arr, delim)
  this_col = arr[col_num]
  num_bits = split(this_col, arr2, "|")
  new_col = this_col
  if (num_bits > 1) {
    new_col = arr2[1]
    for (z=2; z<=num_bits; z++) {
      new_col = new_col"/"arr2[z]
    }
  }
  new_data=new_col
  for (z=2; z<=num_cols; z++) {
    new_data=new_data""delim""arr[z]
  }
  return new_data
}

function fix_NR_NV_col(data, col_num, delim) {
  num_cols = split(data, arr, delim)
  this_col = arr[col_num]
  num_bits = split(this_col, arr2, ",")
  if (num_bits > 1) {
    highest_num = arr2[1]
    for (z=2; z<=num_bits; z++) {
      highest_num = max(highest_num, arr2[z])
    }
    this_col = highest_num
  }
  new_data=arr[1]
  for (z=2; z<=num_cols; z++) {
    if (z==col_num) {
      new_data=new_data""delim""this_col
    } else {
      new_data=new_data""delim""arr[z]
    }
  }
  return new_data
}

function max(a, b) {
    if (a>b) {
        return a
    } else {
        return b
    } fi
}

function abs(a) {
    abs_a = a
    if (a < 0) {
        abs_a = a * -1
    } fi
    return abs_a
}

function avg_of_comma_separated(array_string) {
    avg_of_arr = 0
    num_in_arr = split(array_string,arr,",")
    for (i = 1; i <= num_in_arr; i++) {
      avg_of_arr = avg_of_arr + arr[i]
    }
    avg_of_arr = int(avg_of_arr / num_in_arr)
    return avg_of_arr
}

BEGIN {
    FS = "\t"
    OFS = "\t"
}
{
  chr1 = substr($0, 1, 1)
  if (chr1 == "#") {
    chr6 = substr($0, 1, 6)
    if (chr6 == "#CHROM") {
      print "##INFO=<ID=JOINTCALL_AF,Number=A,Type=Integer,Description=\"Joint-call allele frequency. Is either GATK AF or Platypus FR.\">"
    }
    print $0
  } else {

    # FORMAT field is field 9. See if it has NR or NV.
    # They have comma separated values that must be converted to a single value each.
    # We will not take the average value here. Instead, further below, we will take the maximum value.

    format_nr = -1
    format_nv = -1
    jointcall_af = "."

    num_in_format = split($9, format_array, ":")
    for (j = 1; j <= num_in_format; j++) {
      val = format_array[j]
      if (format_array[j] == "NR") {
        format_nr = j
      }
      if (format_array[j] == "NV") {
        format_nv = j
      }
    }

    # INFO field is field 8. See if it has GATK AF or Platypus FR for joint-call frequency.

    num_in_info = split($8, info_array, ";")
    for (j = 1; j <= num_in_info; j++) {
      val = info_array[j]
      bits = split(val, bits_array, "=")
      if (bits_array[1] == "FR") {
        jointcall_af = bits_array[2]
      }
      if (bits_array[1] == "AF") {
        jointcall_af = bits_array[2]
      }
    }

    $8 = $8";JOINTCALL_AF="jointcall_af

    #if ((format_nr > 0) || (format_nv > 0)) {
    #  for (sample_upto = 10; sample_upto <= NF; sample_upto++) {
    #    new_sample = ""
    #    num_in_sample = split($sample_upto, sample_array, ":")
    #    for (j = 1; j <= num_in_format; j++) {
    #      this_value = sample_array[j]
    #      if ((j == format_nr) || (j == format_nv)) {
    #        this_value = avg_of_comma_separated(this_value)
    #      }
    #      if (new_sample == "") {
    #        new_sample = this_value
    #      } else {
    #        new_sample = new_sample":"this_value
    #      }
    #    }
    #    $sample_upto = new_sample
    #  }
    #}

    # GT needs to have slash (eg. 0/1) instead of bar (eg. 0|1). Fix it if it has bar by changing it to slash.
    # FORMAT field is field 9. See if it has NR or NV.
    # If so, then every sample field (which is fields 10 onwards) needs to be adjusted so that NR and NV contain only 1 value.

    GT_col = find_NR_NV_col($9, "GT", ":")
    NR_col = find_NR_NV_col($9, "NR", ":")
    NV_col = find_NR_NV_col($9, "NV", ":")
    for (i = 10; i <= NF; i++) {
      new_sample = fix_GT_col($i, GT_col, ":")
      new_sample = fix_NR_NV_col(new_sample, NR_col, ":")
      new_sample = fix_NR_NV_col(new_sample, NV_col, ":")
      $i = new_sample
    }

    print $0
  }
}


