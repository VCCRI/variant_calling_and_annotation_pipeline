# Author: Michael Troup and Emma Rath
# Date: Sept 2020
#
# awk program to collapse bedtools intersect output
# to one line per record.  Multiple lines for the same record
# only differ by the last column.  New record will have the 
# last column with values concatinated by ", "
# 1) Assumes input tab delimited
# 2) Assumes input is already in the desired sort order
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

# add_str function - add 2nd string to 1st, only if not already present
function add_str(existing_result, new_element, DELIM) {
    final_result = existing_result
    if ((existing_result == "") || (existing_result == ".")) {
        final_result = new_element
    } else {
        if ((new_element == "") || (new_element == ".")) {
            final_result = existing_result
        } else {
            split(existing_result, existing_result_array, DELIM)
            already_present=0
            for(i=1; i in existing_result_array; i++) {
                existing_element=existing_result_array[i]
                if (new_element==existing_element)
                    already_present=1
                fi
            }
            if(already_present==0)
               final_result = existing_result DELIM new_element
            fi
        }
    }
    return final_result
}

BEGIN {
    FS = "\t"
    OFS = "\t"
    # output delimiter to concatinate last column
    DELIM = ", "
}
NR==1 {
    # print the header
    print
    next
}
{
    # store the last field in a variable
    last_col = $NF
    # make an id from all the fields but the last one
    NF -= 1
    id = $0
    if (NR==2) {
        # special case for the first data record
        # - don't want to print it yet
        last_col_str = last_col
    } else {
        # 2nd or subsequent data line
        # - if new record, then print previous
        if (prev_id != id) {
            # print the previous record
            print prev_id, truncate(last_col_str)
            # start a new last_col_str value
            last_col_str = last_col
        } else {
            # continuation of a record
            last_col_str = add_str(last_col_str, last_col, DELIM)
        }
    }
    prev_id = id
}
END {
    if (NR>1) {
      # print the last record
      print prev_id, truncate(last_col_str)
    }
}

