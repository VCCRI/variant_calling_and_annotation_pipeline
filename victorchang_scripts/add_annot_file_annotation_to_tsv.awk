# Appends a bed file as annotation to tab-delimited file
# After last annotation in INFO field
# Author: Emma Rath
# Date: Jul 2021

# awk -f add_annot_file_annotation_to_tsv.awk -v annot_file="$annot_file" $infile > $outfile

BEGIN {
    FS = "\t"
    OFS = "\t"
    # read annot_file into memory
    line_num = 1
    # starting column for annot_file values to include in output as annotation fields
    START_COL_annot_file = 5
    while (getline < annot_file) {
        # annot_file header like:
        # CHROM	POS	REF	ALT	SpliceAIresult	SpliceAIresult_ALLELE	SpliceAIresult_SYMBOL	SpliceAIresult_DS_AG	SpliceAIresult_DS_ALSpliceAIresult_DS_DG	SpliceAIresult_DS_DL	SpliceAIresult_DP_AG	SpliceAIresult_DP_AL	SpliceAIresult_DP_DG	SpliceAIresult_DP_DL
        # store header names
        if (line_num == 1) {
            annot_file_header[1] = $START_COL_annot_file
            idx = 1
            annot_file_hdr_num_cols = idx
            for (i=(START_COL_annot_file+1); i<=NF; ++i) {
                annot_file_header[++idx] = $i
                annot_file_hdr_num_cols = idx
            }
            ++line_num
        } else {
            # data line
            # do not include true duplicates (all fields same) from annot_file 
            # contraints file
            if (!($0 in seen)) {
                annot_file_key = $1 "|" $2 "|" $3 "|" $4
                # If first time annot_file_key seen, then output will be like:
                # col1val   col2val   col3val
                # allow for different transcripts with same annot_file_key
                # annot_file_key_val will be an associative arrary where the key is the 
                # annot_file_key, and the values are an integer indexed array
                # like: annot_file_key_val[id][1] = col1val
                # annot_file_key_val[id][2] = col2val
                # etc
                if (annot_file_key in annot_file_key_val) {
                    # Already seen annot_file key - so need to append each value separately with
                    # a bar/pipe "|" like:
                    # <gene>|<gene>
                    # <val>|<val>
                    idx = 0
                    for (i=START_COL_annot_file; i<=NF; ++i) {
                        ++idx
                        annot_file_key_val[annot_file_key][idx] = annot_file_key_val[annot_file_key][idx] "|" $i
                    }
                } else {
                    # first time annot_file key seen
                    idx = 0
                    for (i=START_COL_annot_file; i<=NF; ++i) {
                        ++idx
                        annot_file_key_val[annot_file_key][idx] = $i
                    }
                }
                seen[$0]++ # keep track for duplicates
            }
        }
    }
}
(NR==1) {
    # print the header & add the new annotation column headers
    new_hdr = $0
    for (idx=1; idx<=annot_file_hdr_num_cols; ++idx) {
        new_hdr = new_hdr"\t"annot_file_header[idx]
    }
    print new_hdr
    next
}
# recall that input is a tab-delimited file
{
    # id = <chr>|<pos>|<alt>|<ref>
    id = $1 "|" $2 "|" $5 "|" $6
    if (id in annot_file_key_val) {
        new_line = $0
        for (idx=1; idx<=annot_file_hdr_num_cols; ++idx) {
            new_val = annot_file_key_val[id][idx]
            new_line = new_line"\t"new_val
        }
        print new_line
    } else {
        new_line = $0
        this_set_of_values = annot_file_key_val[id]
        for (idx=1; i<=annot_file_hdr_num_cols; ++idx) {
            new_line = new_line"\t."
        }
        print new_line
    }
}

