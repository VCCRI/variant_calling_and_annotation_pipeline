# Author:           Michael Troup
# Date created:     5/5/20
#
# This awk script is was designed to be called by the 
# script: manta_2_annotate_multiple_bedtools_hits.sh
# Intelligent paste of 3 files - two bnd files & an input tsv file
# Generated because sometimes, the 3 files are different lengths
#
# 2 file names passed in as arguments (the bnd files)
# 1 file to process - will print all columns in this file
# add last column from matching record from each of 
# the two arugument files
#
#Usage:
#awk -v left_bnd_file="${left_bnd_intersect_collapsed}" \
#    -v right_bnd_file="${right_bnd_intersect_collapsed}" -f "${awk_script_bnd}" \
#    "${outprefix_dgv_collapsed}" > "${outprefix_segdup}"
#
BEGIN {
    FS = "\t"
    OFS = "\t"
    # input files like:
    #chrom	start	end	right_BND_in_segmental_duplication_region	segmental_duplication_region_hit_by_right_BND
    # 4th col contains key built from an input tsv file:
    #$1 ":" $2 "-" $3 ":" $4 ":" $5
    left_bnd["NA"] = "NA"
    right_bnd["NA"] = "NA"
    left_header = "."
    right_header = "."
    ERROR = "ERROR in combine_bnd_files.awk -" 
    this_nr = 0
    # LEFT - read in the bnd file
    while((getline<left_bnd_file)) {
        ++this_nr
        if (this_nr == 1) {
            # last column
            left_header = $NF
        } else {
            # a pre-existing key is in col 4
            key = $4
            # store the value of the last column (the seg region)
            # in an associative array with the key
            # being the value in col 4
            left_bnd[key] = $NF
        }
    }
    if (this_nr == 0) {
        print ERROR "No records read from file: " left_bnd_file >"/dev/stderr"
        exit 1
    }
    this_nr = 0
    # RIGHT - read in the bnd file
    while((getline<right_bnd_file)) {
        ++this_nr
        if (this_nr == 1) {
            right_header = $NF
        } else {
            key = $4
            right_bnd[key] = $NF
        }
    }
    if (this_nr == 0) {
        print ERROR "No records read from file: " left_bnd_file >"/dev/stderr"
        exit 1
    }
}
# processing the input file line by line - add the new fields
# keep the input tsv columns & add the left & right bnd values
NR==1 {
    #header
    print $0, left_header, right_header
    next
}
{
    #data
    key = $1 ":" $2 "-" $3 ":" $4 ":" $5
    left = "."
    right = "."
    # append the bnd values from the left & right bnd files
    # if the key matches this record
    if (key in left_bnd)
        left = left_bnd[key]
    if (key in right_bnd)
        right = right_bnd[key]
    print $0, left, right
}

