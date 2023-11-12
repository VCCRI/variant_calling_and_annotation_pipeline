# Author: Emma Rath
# Date: 2021-06-16
#
# awk -v min_vaf_for_snvs=0.1 -v min_depth_for_snvs=1 -v min_vaf_for_mnvs=0.2 -v min_depth_for_mnvs=10 -f filter_vcf_snv_and_mnv_vars_for_min_depth_and_min_vaf.awk temp_in.vcf > temp_out.vcf
#
# VPOT filters by sample. When VPOT filtering parameters have Coverage and Hete_Balance, then a sample variant must satisfy both parameters to pass filtering.
#
# This program filters variants. It allows a variant through when it satisfies both minimum depth (of coverage) and minumum vaf (variant allele frequency = alt_depth / total_depth).
#
# This program has one set of filters for SNVs where both REF and ALT are 1-bp in length,
# and another set of filters (usually more stringent) for MNVs where the length of either REF or ALT or both is greater than 1-bp.
#
# This program uses GATK AD and DP fields or Platypus NR and NV fields.
# ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
# ##FORMAT=<ID=NR,Number=A,Type=Integer,Description="Number of reads covering variant location in this sample">
# ##FORMAT=<ID=NV,Number=A,Type=Integer,Description="Number of reads containing variant in this sample">

BEGIN {
  FS = "\t"
  OFS = "\t"
}
{
  chr1 = substr($0, 1, 1)
  if (chr1 == "#") {
    print $0
  } else {

    min_vaf = min_vaf_for_snvs
    min_depth = min_depth_for_snvs
    if ((length($4) > 1) || (length($5) > 1)) { # MNVs have different filtering criteria to SNPs
      min_vaf = min_vaf_for_mnvs
      min_depth = min_depth_for_mnvs
    }

    col_AD = -1
    col_DP = -1
    col_NV = -1
    col_NR = -1
    num_in_format = split($9, format_array, ":")
    for (j = 1; j <= num_in_format; j++) {
      val = format_array[j]
      if (format_array[j] == "NR") {
        col_NR = j
      }
      if (format_array[j] == "NV") {
        col_NV = j
      }
      if (format_array[j] == "AD") {
        col_AD = j
      }
      if (format_array[j] == "DP") {
        col_DP = j
      }
    }

    found_a_sample_variant_that_passes_the_filter = 0

    for (i = 10; i <= NF; i++) {
      depth = -1
      vaf = -1
      split($i, format_array, ":")

      if ((col_AD>-1) && (col_DP>=1)) {
        this_depth = format_array[col_DP]
        if ( this_depth ~ /^[0-9]+$/ ) {
          depth = this_depth
        }
        num_AD = split(format_array[col_AD], AD_array, ",")
        max_AD = -1
        for (k = 2; k <= num_AD; k++) { # the first one is REF, start from the second one which is an ALT
          this_AD = AD_array[k]
          if ( this_AD ~ /^[0-9]+$/ ) {
            if (this_AD > max_AD) {
              max_AD = this_AD
            }
          }
        }
        if ((depth > -1) && (depth >= min_depth)) {
          if ((depth > 0) && (max_AD > -1)) {
            vaf = max_AD/depth
            if (vaf >= min_vaf) {
              found_a_sample_variant_that_passes_the_filter = 1
            }
          }
        }

      } else {

        if ((col_NV>-1) && (col_NR>=1)) {
          this_depth = format_array[col_NR]
          if ( this_depth ~ /^[0-9]+$/ ) {
            depth = this_depth
          }
          num_NV = split(format_array[col_NV], NV_array, ",")
          max_NV = -1
          for (k = 1; k <= num_NV; k++) { # the first one is ALT, there may be a second one and more
            this_NV = NV_array[k]
            if ( this_NV ~ /^[0-9]+$/ ) {
              if (this_NV > max_NV) {
                max_NV = this_NV
              }
            }
          }
          if ((depth > -1) && (depth >= min_depth)) {
            if ((depth > 0) && (max_NV > -1)) {
              vaf = max_NV/depth
              if (vaf >= min_vaf) {
                found_a_sample_variant_that_passes_the_filter = 1
              }
            }
          }
        }

      }
    }

    if (found_a_sample_variant_that_passes_the_filter == 1) {
      print $0
    }
  }
}

