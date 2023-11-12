# Author: Emma Rath
# Date: 2020-12-10
#
# awk -f adjust_vpot_score_for_MNVs.awk input.txt > output.txt
#
# The input file is the vpot priority file in VPOL format:
# Ranking Priority_Score  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    GENE_NAME       19F00715        19W001098
#
# Input will be sorted by CHROM and POS
# When variants overlap with subsequent variants,
# such as in the case of MNVs whose singular SNPs will be present as subsequent variant rows,
# then give that group of variants the ranking and score of the variant that has the highest ranking and score.
# Thus, when the file will be sorted by ranking and score,
# those overlapping variants will still be together as subsequent variants,
# even when they originally had different scores.
#
# This program will look ahead for 20 variants.
# Thus a group of more than 20 overlapping variants will not be correctly treated together as one group.
# Actually, it does a look-behind, reading in 20 lines of variants before outputting the first line.
# When reading in a line, we output the line we read in 20 lines ago. It is stored in prev_line[20].

function test_for_overlap( chrom1, pos1, end1, chrom2, pos2, end2 ) {
  do_they_overlap = 0
  if (chrom1 == chrom2) {
    if ((pos1 <= pos2) && (end1 >= pos2)) { do_they_overlap = 1 }
    if ((pos1 <= pos2) && (end1 >= end2)) { do_they_overlap = 1 }
    if ((pos1 <= end2) && (end1 >= end2)) { do_they_overlap = 1 }
    if ((pos1 >= pos2) && (end1 <= end2)) { do_they_overlap = 1 }
  }
  return do_they_overlap
}

function get_highest_number( string_of_numbers ) {
  num_in_arr = split(string_of_numbers, arr, ",")
  highest_number = arr[1]
  if (num_in_arr > 1) {
    for (x=2; x<=num_in_arr; ++x) {
      if (arr[x] > highest_number) { highest_number = arr[x] }
    }
  }
  return highest_number
}

function set_ranking_and_score( new_ranking, new_score, old_line ) {
  num_in_arr = split(old_line, arr, "\t")
  new_line = new_ranking"\t"new_score
  for (x=3; x<=num_in_arr; ++x) {
    new_line = new_line"\t"arr[x]
  }
  return new_line
}

BEGIN {
  FS = "\t"
  OFS = "\t"
  num_lookahead = 20
  for (i=1; i<=num_lookahead; ++i) {
    prev_line[i] = ""
  }
}

{
  # Output the 20th look-behind row.
  if (prev_line[num_lookahead] != "") {
    print prev_line[num_lookahead]
  }

  # Shift the look-ahead rows by one row,
  # so we will get a new 20th look-behind row which was previously the 19th look-behind row.
  # The row that has just been read in by awk will become the 1st look-behind row.
  for (i=num_lookahead; i>=1; --i) {
    j = i - 1
    prev_line[i] = prev_line[j]
  }
  prev_line[1] = $0

  # From here on we refer to the 20th look-behind row.
  # However, we may not yet have seen 20 rows, and maybe never will if the file is smaller than 20 rows.
  # So figure out what the oldest row is, which will be the row that subsequent rows will be compared to.
  row_to_treat = -1
  for (i=1; i<=num_lookahead; ++i) {
    if (prev_line[i] != "") {
      row_to_treat = i
    }
  }

  # What are the coordinates of the 20th look-behind row,
  # to be compared with other variants to see whether the other variants overlap the 20th look-behind row.
  num_in_arr = split(prev_line[row_to_treat], arr, "\t")
  coords_chrom = arr[3]
  coords_pos = arr[4]
  coords_end = coords_pos + length(arr[6]) - 1
  coords_ranking = arr[1]
  coords_score = arr[2]
  string_of_ranking = coords_ranking
  string_of_score = coords_score

  # Mark which of the 19 look-behind rows overlap with the 20th look-behind row.
  # Of the overlapping rows, what is the highest ranking and score.
  for (i=1; i<=num_lookahead; ++i) {
    overlaps[i] = 0
    if (prev_line[i] != "") {
      num_in_arr = split(prev_line[i], arr, "\t")
      test_chrom = arr[3]
      test_pos = arr[4]
      test_end = test_pos + length(arr[6]) - 1
      test_ranking = arr[1]
      test_score = arr[2]
      overlaps_result = test_for_overlap( coords_chrom, coords_pos, coords_end, test_chrom, test_pos, test_end )
      if (overlaps_result == 1) {
        overlaps[i] = 1
        string_of_ranking = string_of_ranking","test_ranking
        string_of_score = string_of_score","test_score
      }
    }
  }

  # Work out the highest ranking and score.
  # Then set the overlapping variants to have that highest ranking and score.
  highest_ranking = get_highest_number( string_of_ranking )
  highest_score = get_highest_number( string_of_score )
  for (i=1; i<=num_lookahead; ++i) {
    if (overlaps[i] == 1) {
      prev_line[i] = set_ranking_and_score( highest_ranking, highest_score, prev_line[i] )
    }
  }
}

END {
  # Print out the remaining rows.
  # They already have the correct highest ranking and score assigned to them.
  for (i=num_lookahead; i>=1; --i) {
    if (prev_line[i] != "") { print prev_line[i] }
  }
}

