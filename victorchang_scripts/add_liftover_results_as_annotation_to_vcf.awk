# Appends hg19 liftover annotation to vcf file
# After last annotation in INFO field
# Successful liftover file like:
#chrY   2649801 2649802 chrY|2781761|CA|C
#chrY   2649801 2649802 chrY|2781761|CAA|C
# Author: Michael Troup
# Date: Dec 2020
# Original script name: vcf_anno_3.1_hg19_liftover.awk

# awk -f $sw/victorchang_scripts/add_liftover_results_as_annotation_to_vcf.awk -v to_genome="$to_genome" -v liftover_results_bed="$tmp_hg19_out_bed" $infile > $tmp_hg19_out_vcf

# function check for file existance
function check_file(f) {
    # NOTE that system will return 0 if successful ; non-zero on failure
    if (system("[ -f " f " ]")) {
        print "ERROR. TERMINATING.  Invalid file: ", f > "/dev/stderr"
        exit
    }
}

BEGIN {
    FS = "\t"
    OFS = "\t"
    check_file(liftover_results_bed)
    while (getline < liftover_results_bed) {
        if (NF != 4) {
            print "Error in awk liftover script: invalid file format" > "/dev/stderr"
            exit
        }
        # 4th col = id, 3rd col is 1-based start pos (because every entry in the 
        # bed file is for a single base
        # Output like "chr1:1234"
        hg19[$4] = $1 ":" $3
    }
    if (to_genome == ""){
      to_genome = "hg19"
    }
}
/^#/ {
    # print the header & add the new INFO fields
    print
    if ($0 ~/^##INFO=<ID=ALLELE_END/) {
        # add new INFO fields after annovar
        print "##INFO=<ID="to_genome",Number=.,Type=String,Description=\""to_genome" annotation - liftover POS to "to_genome" - added by VCCRI\">"
    }
    next
}
# recall that input is a VCF file
{
    hg19_anno = "."
    # id = <chr>|<pos>|<alt>|<ref>
    id = $1 "|" $2 "|" $4 "|" $5
    if (id in hg19) {
        hg19_anno = hg19[id]
    }
    # add to INFO annotation (col 8 of VCF)
    # - this new annotation will be at the end of the field
    $8 = $8 ";"to_genome"=" hg19_anno
    print
}

