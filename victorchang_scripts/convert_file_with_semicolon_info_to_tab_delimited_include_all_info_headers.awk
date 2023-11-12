# Author: Emma Rath
# Date: 2020-12-10
#
# gawk -f convert_file_with_semicolon_info_to_tab_delimited.awk input.txt > output.txt
#
# The input file has an INFO column with semi-colon-delimited fields, where field is KEY=VALUE
# eg.: BRF=0.56;FR=0.0179;HP=3;HapScore=3;MGOF=228;MMLQ=30;MQ=43.39;NF=161;NR=61;PP=119;QD=0.866972;SC=ATTCCATCCCCTTCCATTCCA;SbPval=0.99
# The rest of the columns are tab-delimited. It may be a VCF file or VPOL file output from vpot priority.
# The header may be just one line or may be the vcf headers of multiple lines starting with a hash.
#
# This program converts the INFO column to tab-delimited values. The column headings are the key.
#
# To obtain all the INFO column headers, even if they don't apear in the first data line,
# this program first reads the entire file to get all the INFO field headers.
# Then this program reads the file again to do that actual processing of each line.

BEGIN {
  FS = "\t"
  OFS = "\t"
  in_hdr = 1
  wrote_header = 0
  prev_line = ""

  info_col = -1
  info_hdrs_string = ""
  idx = 0
  l=0
  while (getline < ARGV[1]) {
    l++
    if (l == 1) {
      for (i=1; i<=NF; ++i) {
        if ($i == "INFO") {
          info_col = i
        }
      }
    } else {
      if (info_col > -1) {
        num1 = split($info_col, arr, ";")
        for (j=1; j<=num1; ++j) {
          split(arr[j], bits, "=")
          new_hdr = bits[1]
          if (info_hdrs_string == "") {
            idx++
            info_hdrs[idx] = new_hdr
            info_hdrs_string = "got_first_value"
          } else {
            already_got_hdr = 0
            for (k=1; k<=idx; ++k) {
              if (info_hdrs[k] == new_hdr) {
                already_got_hdr = 1
              }
            }
            if (already_got_hdr == 0) {
              idx++
              info_hdrs[idx] = new_hdr
            }
          }
        }
      }
    }
  }
  close(ARGV[1])
  # sort fields so they always appear in same order regardless of what info fields are not present in first data row
  asort(info_hdrs)
  # print out the info headers for debuggin
  #for (i=1; i<=idx; ++i) { print i, info_hdrs[i] }
  # make an associative array, and make an initialised associative_array
  for (i=1; i<=idx; ++i) {
    val = info_hdrs[i]
    info_hdrs_assoc[val] = i
    info_hdrs_init[val] = "."
  }
}

{
  if (in_hdr == 1) { # We are still in header rows, unless this is the first data row.

    if (NR > 1) {
      chr1 = substr($1, 1, 1)
      if (chr1 != "#") {
        in_hdr = 0
      }
    }
    if (in_hdr == 1) {
      prev_line = $0
    }
  }

  if (in_hdr == 0) { # We are now in the data rows.

    # We have most of the header fields saved in prev_line, and we already got the INFO header fields from previously running through the data.
    if (wrote_header == 0) {

      split(prev_line, prev, "\t")
      for (i=1; i<=NF; ++i) {

        if (i>1) { printf OFS } # Put a tab before every column except for the first column

        if (i == info_col) {
          # This is the INFO field, so write out its keys fields as tab-delimited headers
          for (j=1; j<=idx; ++j) {
            this_hdr = info_hdrs[j]
            if (j>1) { printf OFS } # Put a tab before every column except for the first column
            printf this_hdr
          }

        } else {
          # This is not the INFO field, so just write out this header
          printf prev[i]
        }
      }
      printf "\n"
      wrote_header = 1
    }

    # Now write out this data record, splitting up its INFO fields
    for (i=1; i<=NF; ++i) {

      if (i>1) { printf OFS } # Put a tab before every column except for the first column

      if (i == info_col) {
          # This is the INFO field, so first get all the subfield values
          # Then write out all the values as tab-delimited
          for (j=1; j<=idx; ++j) {
             this_key = info_hdrs[j]
             this_data_init[this_key] = info_hdrs_init[this_key]
          }
          num = split($i, arr, ";")
          for (j=1; j<=num; ++j) {
            num2 = split(arr[j], bits, "=")
            this_key = bits[1]
            new_data = "."
            if (num2 > 1) { new_data = bits[2] }
            this_data_init[this_key] = new_data
          }
          # Then write out all the values as tab-delimited
          for (j=1; j<=idx; ++j) {
            this_key = info_hdrs[j]
            this_data = this_data_init[this_key]
            if (j>1) { printf OFS } # Put a tab before every column except for the first column
            {printf "%s", this_data}
          }

     } else {
        # This is not the INFO field, so just write out this data
        printf $i
      }
    }
    printf "\n"
  }
}

