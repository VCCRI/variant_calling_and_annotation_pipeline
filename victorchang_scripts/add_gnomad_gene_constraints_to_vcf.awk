# Appends gnomad gene constraints annotation to vcf file
# After last annotation in INFO field
# Input variable: bed - file with intersections of VCF site & gnomad constraints site
# Input variable: gnom - gnomad constraints file
# Input variable: chr - in form "chr1", "chrY", etc
# Input FILE: VCF file
# Author: Michael Troup
# Date: Dec 2020
# Original script name: vcf_anno_3.2_gnomad_gene_constraints.awk

# function check for file existance
function check_file(f) {
    # NOTE that system will return 0 if successful ; non-zero on failure
    if (system("[ -f " f " ]")) {
        print "ERROR. TERMINATING.  Invalid file: ", f > "/dev/stderr"
        exit 1
    }
}

BEGIN {
    FS = "\t"
    OFS = "\t"
    if (length(chr) < 1) {
        print "Error in awk gene constraints script: missing variable - chr" > "/dev/stderr"
        exit 1
    }
    check_file(bed)
    # bed input file is from positive matches of VCF sites with gnomad constraints
    # file.  Input looks like:
    #chrY	12538348	12538348	chrY|12538348|C|T	chrY	12537649	12860839	chrY|12537649|12860839|USP9Y	0
    #chrY	12538979	12538979	chrY|12538979|T|A	chrY	12537649	12860839	chrY|12537649|12860839|USP9Y	0
    # Col 4 is VCF id ; col 8 is gnomad constraints id
    while (getline < bed) {
        if (NF < 8) {
            print "Error in awk gene constraints script: invalid file format" > "/dev/stderr"
            exit 1
        }
        # could be multiple lines with same VCF id (col 4), due to different transcripts
        # (only catching diff transcript on diff gene here)
        if ($4 in vcf_gnomad)
            vcf_gnomad[$4] = vcf_gnomad[$4] "," $8
        else
            vcf_gnomad[$4] = $8
    }
    # read gnomad constraints into memory
    check_file(gnom)
    line_num = 1
    # starting column for gnomad constraints values to include in INFO fields
    START_COL_GNOM = 4
    while (getline < gnom) {
        # gnomad constraints header like:
        # chrom	start	end	gnomad_constraints_gene	gnomad_constraints_syn_z	gnomad_constraints_oe_syn	gnomad_constraints_oe_syn_lower	gnomad_constraints_oe_syn_upper	gnomad_constraints_mis_z	gnomad_constraints_oe_mis	gnomad_constraints_oe_mis_lower	gnomad_constraints_oe_mis_upper	gnomad_constraints_pLI	gnomad_constraints_oe_lof	gnomad_constraints_oe_lof_lower	gnomad_constraints_oe_lof_upper
        # store header names + create a default entry for "no match"
        if (line_num == 1) {
            gnom_header[1] = $4
            default_info[1] = gnom_header[1] "=."
            idx = 1
            for (i=5; i<=NF; ++i) {
                gnom_header[++idx] = $i
                default_info[idx] = $i "=."
            }
            ++line_num
        } else {
            # data line
            # only need to read in for this chr
            if ($1 == chr) {
                # do not include true duplicates (all fields same) from gnomad 
                # contraints file
                if (!($0 in seen)) {
                    gnom_key = $1 "|" $2 "|" $3 "|" $4
                    # If first time gnom_key seen, then output will be like:
                    # gnomad_constraints_gene=<gene>;gnomad_constraints_syn_z=<val>;...
                    # allow for different transcripts with same gnom_key
                    # gnom_key_val will be an associative arrary where the key is the 
                    # gnomad constraint id, and the values are an integer indexed array
                    # like: gnom_key_val[id][1] = gnomad_constraints_gene=<gene>
                    # gnom_key_val[id][2] = gnomad_constraints_syn_z=<val>
                    # etc
                    if (gnom_key in gnom_key_val) {
                        # Already seen gnomad key - so need to append each value separately with
                        # a comma "," like:
                        # <gene>,<gene>
                        # <val>,<val>;...
                        idx = 0
                        for (i=START_COL_GNOM; i<=NF; ++i) {
                            ++idx
                            gnom_key_val[gnom_key][idx] = gnom_key_val[gnom_key][idx] "," $i
                        }
                    } else {
                        # first time gnomad key seen
                        # Don't want to include name of info field yet
                        idx = 0
                        for (i=START_COL_GNOM; i<=NF; ++i) {
                            ++idx
                            gnom_key_val[gnom_key][idx] = $i
                        }
                    }
                    seen[$0]++ # keep track for duplicates
                }
            }
        }
    }
}
/^#/ {
    # print the header & add the new INFO fields 
    if ($0 ~/^#CHROM/) {
        # add new INFO fields after annovar
        for (i=1; i<=length(gnom_header); ++i) {
            print "##INFO=<ID=" gnom_header[i] ",Number=.,Type=String,Description=\"" gnom_header[i] " annotation - gnomAD gene constraints - added by VCCRI\">"
        }
    }
    print
    next
}
# recall that input is a VCF file
{
    # id = <chr>|<pos>|<alt>|<ref>
    id = $1 "|" $2 "|" $4 "|" $5
    if (id in vcf_gnomad) {
        #print "id=" id > "/dev/stderr" #TODO debugging
        # check if multiple gnomad values for the one vcf id
        r = split(vcf_gnomad[id], gnomad_ids, ",")
        if (!(gnomad_ids[1] in gnom_key_val)) {
            print "Error in awk gene constraints script: invalid gnomad key: " gnomad_ids[1] > "/dev/stderr"
            print "id=" id > "/dev/stderr"
            print "i=1" > "/dev/stderr"
            print "vcf_gnomad[id]=<" vcf_gnomad[id] ">" > "/dev/stderr"
            exit 1
        }
        # TODO - check length of arrays are as expected
        # process first gnomad id - prefix with "<header>="
        for (i=1; i<=length(gnom_key_val[gnomad_ids[1]]); ++i) {
            extra_anno[i] = gnom_header[i] "=" gnom_key_val[gnomad_ids[1]][i]
        }
        # process second & subsequent gnomad ids (if any) - do NOT prefix with <header>=
        for (i=2; i<=length(gnomad_ids); ++i) {
            if (!(gnomad_ids[i] in gnom_key_val)) {
                print "Error in awk gene constraints script: invalid gnomad key: " gnomad_ids[i] > "/dev/stderr"
                print "id=<" id ">"> "/dev/stderr"
                print "i=" i > "/dev/stderr"
                print "vcf_gnomad[id]=" vcf_gnomad[id] > "/dev/stderr"
                exit 1
            }
            for (j=1; j<=length(gnom_key_val[gnomad_ids[i]]); ++j) {
                extra_anno[j] = extra_anno[j] "," gnom_key_val[gnomad_ids[i]][j]
            }
        }
    } else {
        for (i=1; i<=length(default_info); ++i)
            extra_anno[i] = default_info[i]
    }
    # add to INFO annotation (col 8 of VCF)
    # - this new annotation will be at the end of the field
    for (i=1; i<=length(extra_anno); ++i)
        $8 = $8 ";" extra_anno[i]
    print
}

