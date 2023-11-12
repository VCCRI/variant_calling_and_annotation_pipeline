# Author: Emma Rath
# Modified version of script by Michael Troup
# Date: 2020-07-17
#
# awk -v key -v collapse -f collapse_manta_bedtools_SVs_to_unique_by_column.awk my_input.tsv > my_output.tsv
#
# awk -v key=1,2,3,4,5 -v collapse=12,13,14,15,16,17,18,19,20,21,22,23,24 -f collapse_manta_bedtools_SVs_to_unique_by_column.awk test_collapse1.txt > test_collapse1_collapsed.txt
#
# awk program to collapse bedtools intersect output
# to one line per record.  Multiple lines for the same record
# only differ by the last one or more columns.  New record will have the 
# last columns with values concatenated by ", "
# 1) Assumes input tab delimited
# 2) Assumes input is already in the desired sort order
#
# key is the comma-separated list of the columns that make up the key
# collapse is the comma-separated list of the columns that need to be collapsed to one row per key
# They must be the last one or more columns.
#
###########################################################
# 3) NOTE *** Assumes input file has header line *** NOTE #
###########################################################
# 
# Example
# Input file:
#one two three
#1 2 3
#1 2 4
#1 2 5
#1 3 1
#
# Output:
#one two three
#1 2 3, 4, 5
#1 3 1
#
# collapse_manta_bedtools_SVs_to_unique_row_and_value.awk checks for duplicate values and takes too long to run
# collapse_manta_bedtools_SVs_to_unique.awk does not check for duplicate values and for this field we do not need to check
# collapse_manta_bedtools_SVs_to_unique_by_column.awk does not check for duplicate values and collapses multiple last columns, not just the one last column



# truncate function - returns string truncated to certain length
function truncate(str) {
  MAX_STR = 9999
  if (length(str) > MAX_STR)
    return substr(str,1,MAX_STR) " plus more"
  else
    return str
}
function array_length(a, i, k) {
  k = 0
  for (i in a) k++
  return k
}
BEGIN {
  FS = "\t"
  OFS = "\t"
  # output delimiter to concatenate values in columns that need to be collapsed
  DELIM = ", "
  # initialise work variables
  split(key, id_array, ",")
  split(collapse, collapse_array, ",")
  for (i = 1; i <= array_length(collapse_array, 0, 0); i++) {
   output_collapse_string[i] = ""
  }
  prev_id = ""
  prev_line = ""
}
NR==1 {
  # print the header
  print
  next
}
{
  # make an id from all the fields of the key
  key_col = id_array[1]
  id = $key_col
  if (array_length(id_array, 0, 0) > 1) {
    for (i = 2; i <= array_length(id_array, 0, 0); i++) {
      key_col = id_array[i]
      id = id "," $key_col
    }
  }

  # if this is a new id, then print the one collapsed record for the previous id
  if (prev_id != id) {
    if (prev_id != "") {
      # build the previous record. the last one or more column are the running concatenations of those columns.
      split(prev_line, prev_line_array, FS)
      new_prev_line = prev_line_array[1]
      j = 0
      for (i = 2; i <= NF; i++) {
        if (i < collapse_array[1]) {
          new_prev_line = new_prev_line OFS prev_line_array[i]
        } else {
          j = j + 1
          new_prev_line = new_prev_line OFS output_collapse_string[j]
        }
      }
      # print the previous record
      print new_prev_line
      # start a new running concatenation of the last column(s)
      for (i = 1; i <= array_length(collapse_array, 0, 0); i++) {
        output_collapse_string[i] = ""
      }
    }
  }

  # do the concatenations of the fields that need to be concatenated
  for (i = 1; i <= array_length(collapse_array, 0, 0); i++) {
    this_collapse_col = collapse_array[i]
    if ((output_collapse_string[i] == ".") || (output_collapse_string[i] == "")) {
      output_collapse_string[i] = $this_collapse_col
    } else {
      if (($this_collapse_col == ".") || ($this_collapse_col == "")) {
        output_collapse_string[i] = output_collapse_string[i]
      } else {
        output_collapse_string[i] = output_collapse_string[i] DELIM $this_collapse_col
      }
    }
  }
  prev_id = id
  prev_line = $0
}
END {
  # print the last record
  if (prev_id != "") {
    # build the previous record. the last one or more column are the running concatenations of those columns.
    new_prev_line = $1
    j = 0
    for (i = 2; i <= NF; i++) {
      if (i < collapse_array[1]) {
        new_prev_line = new_prev_line OFS $i
      } else {
        j = j + 1
        new_prev_line = new_prev_line OFS output_collapse_string[j]
      }
    }
    if (NR>1) {
      # print the previous record
      print new_prev_line
    }
  }
}


