# Appends an annotation text to a vpot variant in a vpot tab-delimited file when the bed file overlaps with with the variant
# Author: Emma Rath
# Date: Jul 2021

# awk -f add_annot_region_file_annotation_to_tsv.awk -v annot_file="$annot_file" -v annot="BigHomozygRun" $infile > $outfile

# infile is in tab-delimited format, which may or may not have the END column. Use END if it is present, else use POS+length(REF):
#
# #CHROM  POS     END     ID      REF     ALT     QUAL
# chr17   14236571        14236574        .       CCAA    TCAC
# chr1    3186271 3186271 .       C       T       558.11  PASS
# chr7    37857299        37857299        .       T       C
#
# or
#
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER >
# chr19   5711980 .       C       T       755.19  PASS    AC=1;AF>
# chr1    237890634       .       A       G       315.19  PASS   >
# chr2    158815810       rs138143887     C       T       1098.44>


# annot_file of regions that might overlap vpot variants is in bed format:
#
# homozygosity_run_chrom	hrun_start	hrun_end	hrun_length	hrun1	hrun2	hrun3	hrun4
# chr1	858567	869379	10812	0	22	0	NA
# chr1	1037956	1047614	9658	0	28	0	NA
# chr1	1318756	1361438	42682	0	68	4	NA
#
# or
#
# chrom	start	end
# chr6	6382383	11905641
# chr8	124990054	135826193
# chr9	83224216	92173340


BEGIN {
    FS = "\t"
    OFS = "\t"
    # read annot_file into memory
    num_annot_lines = 0
    seen_hdr = 0
    while (getline < annot_file) {
        # annot_file header like one of the following:
        # homozygosity_run_chrom	hrun_start	hrun_end	hrun_length	hrun1	hrun2	hrun3	hrun4
        # chrom	start	end
        if (seen_hdr == 1) { # skip first line, and store second line as first line
            num_annot_lines = num_annot_lines + 1
            annot_chrom[num_annot_lines] = $1
            annot_start[num_annot_lines] = $2
            annot_end[num_annot_lines] = $3
        }
        seen_hdr = 1
    }
}
(NR==1) {
    # print the header & add the new column header
    print $0, annot
    # do we have an END field? do we have a REF field?
    end_col=-1
    ref_col=-1
    if (($3 == "END") || ($3 == "End") || ($3 == "end")) {
      end_col = 3
    }
    if (($4 = "REF") || ($4 = "Ref") || ($4 = "ref")) {
      ref_col = 4
    }
    if (($5 = "REF") || ($5 = "Ref") || ($5 = "ref")) {
      ref_col = 5
    }
    chrom_col = 1
    pos_col = 2
    next
}
# recall that input is a tab-delimited file
{
    end_val = $pos_col # end is same as start/pos, unless we have an END field or REF field to give the end value
    if (end_col != -1) {
      end_val = $end_col
    } else {
      if (ref_col != -1) {
        end_val = $pos_col + length($ref_col) - 1
      }
    }
    variant_overlaps_annot = 0
    if (num_annot_lines > 0) {
      look_for_overlaps = 1
      idx = 1
      while (look_for_overlaps == 1) {
        if ($chrom_col == annot_chrom[idx]) {
          if ((annot_start[idx] <= $pos_col) && ($pos_col <= annot_end[idx])) {
            variant_overlaps_annot = 1
          } else if (($pos_col <= annot_start[idx]) && (end_val >= annot_start[idx])) {
            variant_overlaps_annot = 1
          } else if ((end_val >= annot_end[idx]) && ($pos_col <= annot_end[idx])) {
            variant_overlaps_annot = 1
          }
        }
        if (variant_overlaps_annot == 1) {
          look_for_overlaps = 0
        }
        idx = idx + 1
        if (idx > num_annot_lines) {
          look_for_overlaps = 0
        }
      }
    }
    if (variant_overlaps_annot == 1) {
      print $0, annot
    } else {
      print $0, "."
    }
}

