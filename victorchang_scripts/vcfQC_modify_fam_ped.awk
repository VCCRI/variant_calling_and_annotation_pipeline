# Script written by Michael Troup, 2019 May
# This file is called by vcf_qc_2_plink.sh
# awk file to modify a .fam file generated by vcftools
# with values from a ped file
# The reason for this file (& not just replacing with headerless ped file)
# is to retain the same order as the .fam file - in case this matters
# .fam files assumed to have the following format:
#A text file with no header line, and one line per sample with the following six fields:
#Family ID ('FID')
#Individual ID ('IID'; cannot be '0')
#Individual ID of father ('0' if father isn't in dataset)
#Individual ID of mother ('0' if mother isn't in dataset)
#Sex code ('1' = male, '2' = female, '0' = unknown)
#Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
#   - no header line
#   - 6 columns
#   - <family id> <subject id> <parent 1> <parent 2> <gender> <some other column>
# ped file assumed to have the following format
#   - header line
#   - Family    Subject Father  Mother  Sex Phenotype
# Required parameters: fam & ped - names of input files
# New .fam file is printed to stdout - calling program pipes to file
BEGIN {
    # check that the required arguments have been supplied
    if (length(fam) < 1 || length(ped) < 1) {
        print "TERMINATING from modify_fam_ped.awk - missing arugment"
        exit 1
    }
    # check that the files exist
    cmd = "[ -f \"" fam "\" ]"
    if (system(cmd) != 0) {
        print "TERMINATING from modify_fam_ped.awk. Can't find fam file: " fam
        exit 1
    }
    cmd = "[ -f \"" ped "\" ]"
    if (system(cmd) != 0) {
        print "TERMINATING from modify_fam_ped.awk. Can't find ped file: " ped
        exit 1
    }
    # read in the ped file - first line is header
    getline < ped
    while (getline < ped)
        a[$2] = $1 " " $2 " " $3 " " $4 " " $5 " " $6
    close(ped)
    if (length(a) < 1) {
        print "TERMINATING from modify_fam_ped.awk. Empty ped file: " ped
        exit 1
    }
    # print the new .fam output to stdout
    while (getline < fam) {
        if (!($2 in a)) {
            print "TERMINATING from modify_fam_ped.awk. Can't find subject in ped file: " $2
            exit 1
        }
        print a[$2]
    }
}