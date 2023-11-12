# awk -f convert_long_variants_into_multiple_gnomad_MNV_format_variants.awk infile > outfile

# Gnomad provided allele frequencies for multi-nucleotide variants/polymorphisms (MNVs/MNPs) at
# https://gnomad.broadinstitute.org/downloads
#
# There are 2bp and 3bp MNVs:
# head hg38_gnomadMNV.txt
# #Chr	Start	End	Ref	Alt	gnomadMNV_AF_mnv	gnomadMNV_n_indv	gnomadMNV_n_homhom	gnomadMNV_mnv_consequence
# chr1	69511	69513	ACA	GCG	2.8495873797474126e-05	4.0	1	missense_variant
# chr1	69946	69947	GA	AG	4.3301290378453275e-05	1.0	0	missense_variant
# chr1	138593	138595	GGC	TGT	1.0627218431847649e-05	1.0	0	missense_variant
# chr1	138801	138802	CT	TG	7.334819856824316e-06	1.0	0	missense_variant
#
# There are MNVs that are 1bp SNPs with a separation distance of 1 to 10 bp
# (adjacent SNPs are separated by 1bp, there are 9bp inbetween 2 SNPs that have a separation distance of 10bp):
# head hg19_gnomadMNVdist.txt
# #Chr	Start	End	Ref	Alt	gnomadMNVdist_AF_mnv	gnomadMNVdist_n_indv	gnomadMNVdist_n_homhom	gnomadMNVdist_mnv_consequence
# 1	10247	10248	TA	CT	4	0	2.82773e-05	0
# 1	10247	10257	TAAACCCTAAA	CAAACCCTAAC	1	0	7.06934e-06	0
# 1	10248	10257	AAACCCTAAA	TAACCCTAAC	1	0	7.06934e-06	0
# 1	10250	10257	ACCCTAAA	CCCCTAAC	21	2	0.000148456	1.41387e-05
#
# The input to this program is a VCF with MNV variants.
# This program splits up the MNV variants into variants that will match the format of how we have formatted gnomad MNVs for annovar annotation.
# Thus, the split-up MNV variants that are output from this program will be annotated with their gnomad MNV frequency, and no annotation means that there are zero of this MNV in gnomad.
# Thhis program does not output the original long MNP because that is not a format that will match gnomad MNVs.
# The formats output by this program, that match gnomad MNVs, are the following:
#  * 2 SNPs side by side that make a 2bp MNV
#  * 3 SNPs side by side that make a 3bp MNV
#  * 2 SNPs with 1 to 9 bp in-between, such that the in-between bp match the reference sequence (no SNPs in-between). 
#
# When looking for 2bp and 3bp MNVs, this program looks from left to right, and also from right to left in case there are indels and the genes involved would be on the reverse strand.
# 
# The Gnomad distance MNVs frequencies do not make assumptions about whether there are SNPs in-between the 2 SNPs that are distance-bp apart.
# However, the annovar reference table that we created does make an assumption that there are no SNPs in-between, and has REF and ALT sequences according to that assumption.
# Thus, this programs only looks for distance MNVs that have no SNPs in-between (because those are the only ones that will match our annovar reference table).
# This program looks from left to right when looking for distance MNVs (because looking from right to left might create too many misleading MNVs).
#
# This program may end up outputting the same variant more than once from the same input variant.
# After this program is run, it is expected that output will be sorted and unique records keep and duplicate records discarded.
# This program does not make sure that it does NOT output the same variant more than once.
#
# This program may end up outputting the same variant more than once from two different input variants
# that themselves were originally from the same multi-allelic variant that was split into 2 variants before calling this program.
# The rest of the non-sample vcf fields will have the same values.
# The vcf sample fields may have the same or different values.
# The sort and uniq operations carried out on the output from this program will remove any such duplicates when all the vcf and samples fields are the same.
# The sort and uniq operations will not remove these duplicates when the sample fields are not the same.
# (eg. when sample GT in one variant = ./0 and sampt GT in another variant = 1/0)
#
# Please note that platypus may have produced 2 variants in the same region (first var and second var),
# and when this program outputs a subset of first var it may end up being the same as second var.
# In this case, the sort and uniq step after this program will not remove the duplicat variant
# because the other fields in the record (info, qual, etc) will not be identical.
# eg.
# input:
# chr10   10547   .       GCTGATCAGGACG   CCTGATGAGGACA   1319    MQ;alleleBias   BRF=0.51;FR=0.1399
# chr10   10553   .       CAGGACG AAGGACA 		1175	MQ;alleleBias   BRF=0.51;FR=0.1417
# output:
# chr10   10553   .       CAGGACG GAGGACA 		1319    MQ;alleleBias   BRF=0.51;FR=0.1399
# chr10   10553   .       CAGGACG AAGGACA 		1175    MQ;alleleBias   BRF=0.51;FR=0.1417
#
# Here is an example where this program copies an MNP and does not copy a non-MNV.
# eg. output from platypus:
# chr10   10980   .       GGG     CGA,CGG 555     badReads;MQ     BRF=0.74;FR=0.1865
# after vt normalize we have:
# chr10   10980   .       GGG     CGA     555     badReads;MQ     BRF=0.74;FR=0.1865
# chr10   10980   .       GGG     CGG     555     badReads;MQ     BRF=0.74;FR=0.1865
# This program does not output the second var which is really only one SNP and is not a true MNP. It outputs only:
# chr10   10980   .       GGG     CGA     555     badReads;MQ     BRF=0.74;FR=0.1865
#
# Here is an example of this program outputting a subset MNV from a larger input MNV only when the ALT base-pairs in-between match the REF base-pairs.
# When the ALT base-pairs in-between do not match the REF base-pairs, then the MNV is not outputted.
# input:
# chr10   52855   .       AGTCATGTG       CATCATGTT       137     alleleBias      BRF=0.15;FR=0.0551
# chr10   52855   .       AGTCATGTG       CGTCATGTT       137     alleleBias      BRF=0.15;FR=0.0551
# output:
# chr10   52855   .       AG      CA      137     alleleBias      BRF=0.15;FR=0.0551
# chr10   52855   .       AGTCATGTG       CGTCATGTT       137     alleleBias      BRF=0.15;FR=0.0551
# chr10   52856   .       GTCATGTG        ATCATGTT        137     alleleBias      BRF=0.15;FR=0.0551
#
# Gnomad-MNV has produced frequencies for MNVs a certain distance apart (from 1 to 10 bp distance),
# regardless of whether the ALT base-pairs in-between match the REF base-pairs or not.
# We use annovar to annovate MNVs, using an annovar reference table that we created from the gnomad-MNV frequencies.
# Annovar requires the full REF and ALT to be specied.
# For MNVs a certain distance apart, we have added a variant containing reference genome sequence in-between the 2 SNPs, for both REF and ALT.
# Thus, a sample variant processed by this program that contains non-matching base-pairs in-between the 2 SNPs that compose the MNV
# will not match the annovar gnomad-MNV annotation reference table.
# Thus, we do not output such MNVs, because they won't match the gnomad-MNV distance frequencies annotation reference table.
#
# This program will output subsetted MNVs that are 2bp or 3bp long, where the 1st and last bp of REF are both different to ALT,
# regardless of whether the ALT base-pairs on either side of the 2bp or 3bp MNV are different to REF.
# This is to ensure that any 2bp or 3bp MNVs that would be wholly inside a codon will not be missed and will be outputted,
# in case they result in a stop-codon, to make that variant stop-codon apparent.
# For the subsetted 2bp MNVs, perhaps we could change this to output it only if the base-pairs on either side are the same in REF and ALT,
# because when they are not the same, they will be output as subsetted 3bp MNV(s).
# Does gnomad-MNV frequencies put an entry in the 2bp file when it is part of a 3bp MNV and will be in the 3bp MNV file?
# This should be verified. If the gnomad-MNV frequency 2bp file contains 2bp MNVs only when there are not SNPs in either side,
# then perhaps we should do the same with output subsetted 2bp MNVs, so that it picks up the gnomad-MNV frequency annotation in the same situation only.
# However, the gnomad-MNV Wang et al. 2020 paper says that it called MNVs in the following way:
# "For MNV discovery, we exhaustively looked for variants that appear in the same individual, in cis,
#  and within 2 bp distance for the exome exome data set and 10 bp distance for the genome data set,
#  using the hail window_by_locus function."
# At https://gnomad.broadinstitute.org/downloads, the readme.md for the "Multi-nucleotide variants (MNVs)" says the following:
# "For example, when
#  locus=22:16056126
#  refs="G,C"
#  alts="C,A" and
#  d=2,
#  the MNV is GcC to CcA, when there is no SNP in the position 22:16056127.
#  (c is the reference allele at position 22:16056127).
#  (Caution: Since we did not explicitly call the MNVs consisting of more than 3 SNPs other than the coding region,
#  we do not exclude the possibility that the MNV is in reality part of larger substitution.
#  For example above, we do not exclude the possibility that the MNV is actually part of GcC to CTA by this file alone.
#  Users can visit the gnomAD browser for manual check.)"
# Thus, gnomad-MNVs that fall outside of coding region may have other SNPs inside the bookend SNPs of the MNV.
# When the gnomad source data contained a string of more than 3-bp where every ALT bp is different to its REF bp,
# has gnomad produced the frequencies of the 2-bp and 3-bp subsetted MNVs, even though there are SNPs on either side?
# We probably need to study the hail window_by_locus function to find out.
# The 3bp gnomad-MNV frequencies do contain 3bp MNVs where all 3 ALT base-pairs are different to their corresponding bp in the REF.


function min(n1, n2) {rtn=n1; if (n2<n1) rtn=n2; return rtn}

BEGIN {
  FS="\t";OFS="\t";
}
{
  char1 = substr($0,1,1)
  if (char1 != "#") {

    # look for 2bp MNVs going from left to right, that are only 2bp long (the nucleotides on either side are REF)
    # this is also an MNV made of 2 SNPs that are a distance of 1bp apart
    # this will include variants whose REF and ALT are originally only 2bp long, and sub-variants in an originally longer variant
    # this 2bp MNV may be embedded in a platypus variant that is an indel
    if (length($4) >= 2) {
      num_char=min(length($4), length($5))
      for (i=2; i<=num_char; i++) {
        r1="N"; if (i>2) { r1=substr($4,(i-2),1) }
        r2=substr($4,(i-1),1)
        r3=substr($4,i,1)
        r4="N"; if (i<num_char) { r4=substr($4,(i+1),1) }
        a1="N"; if (i>2) { a1=substr($5,(i-2),1) }
        a2=substr($5,(i-1),1)
        a3=substr($5,i,1)
        a4="N"; if (i<num_char) { a4=substr($5,(i+1),1) }
        if ((r1==a1)&&(r2!=a2)&&(r3!=a3)&&(r4==a4)) {
          printf $1 OFS ($2+i-2) OFS $3 OFS r2 r3 OFS a2 a3; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
        }
      }
    }

    # look for 2bp MNVs going from right to left, that are only 2bp long (the nucleotides on either side are REF)
    # this 2bp MNV may be embedded in a platypus variant that is an indel
    if (length($4) >= 2) {
      if (length($4) != length($5)) {
        addto4=0; if (length($4) > length($5)) { addto4=length($4)-length($5) }
        addto5=0; if (length($5) > length($4)) { addto5=length($5)-length($4) }
        if (num_char>=3) {
          for (i=num_char; i>=2; i--) {
            r1="."; if (i>2) { r1=substr($4,(i+addto4-2),1) }
            r2=substr($4,(i+addto4-1),1)
            r3=substr($4,i+addto4,1)
            r4="."; if (i<num_char) { r4=substr($4,(i+addto4+1),1) }
            a1="."; if (i>2) { a1=substr($5,(i+addto5-2),1) }
            a2=substr($5,(i+addto5-1),1)
            a3=substr($5,i+addto5,1)
            a4="."; if (i<num_char) { a4=substr($5,(i+addto5+1),1) }
            if ((r1==a1)&&(r2!=a2)&&(r3!=a3)&&(r4==a4)) { # 2 adjacent SNPs makes a 2bp MNV that has 2 SNPs in it
              printf $1 OFS ($2+i-2+addto4) OFS $3 OFS r2 r3 OFS a2 a3; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
            }
          }
        }
      }
    }

    # look for 2bp MNVs regardless of whether nucleotides match the REF on either side
    # only look in platypus variants having same length REF and ALT, so that we do not get a lot of spurious sub-variants from stuttering
    if ((length($4) >= 2) && (length($4) == length($5))) {
      num_char=length($4)
      for (i=2; i<=num_char; i++) {
        r1=substr($4,(i-1),1)
        r2=substr($4,i,1)
        a1=substr($5,(i-1),1)
        a2=substr($5,i,1)
        if ((r1!=a1)&&(r2!=a2)) {
          printf $1 OFS ($2+i-2) OFS $3 OFS r1 r2 OFS a1 a2; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
        }
      }
    }

    # look for 3bp MNVs going from left to right, that are only 3bp long (the nucleotides on either side are REF)
    # this will include variants whose REF and ALT are originally only 3bp long, and sub-variants in an originally longer variant
    # this 3bp MNV may be embedded in a platypus variant that is an indel
    if (length($4) >= 3) {
      for (i=3; i<=num_char; i++) {
        r1="."; if (i>3) { r1=substr($4,(i-3),1) }
        r2=substr($4,(i-2),1)
        r3=substr($4,(i-1),1)
        r4=substr($4,i,1)
        r5="."; if (i<num_char) { r5=substr($4,(i+1),1) }
        a1="."; if (i>3) { a1=substr($5,(i-3),1) }
        a2=substr($5,(i-2),1)
        a3=substr($5,(i-1),1)
        a4=substr($5,i,1)
        a5="."; if (i<num_char) { a5=substr($5,(i+1),1) }
        if ((r1==a1)&&(r2!=a2)&&(r3!=a3)&&(r4!=a4)&&(r5==a5)) { # 3 adjacent SNPs makes a 3bp MNV that has 3 SNPs in it
          printf $1 OFS ($2+i-3) OFS $3 OFS r2 r3 r4 OFS a2 a3 a4; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
        }
        if ((r1==a1)&&(r2!=a2)&&(r3==a3)&&(r4!=a4)&&(r5==a5)) { # SNP then non-SNP then SNP makes a 3bp MNV that has 2 SNPs in it not 3 SNPs
          printf $1 OFS ($2+i-3) OFS $3 OFS r2 r3 r4 OFS a2 a3 a4; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
        }
      }
    }

    # look for 3bp MNVs going from right to left, that are only 3bp long (the nucleotides on either side are REF)
    # this 3bp MNV may be embedded in a platypus variant that is an indel
    if (length($4) >= 3) {
      if (length($4) != length($5)) {
        if (num_char>=4) {
          for (i=num_char; i>=3; i--) {
            r1="."; if (i>3) { r1=substr($4,(i+addto4-3),1) }
            r2=substr($4,(i+addto4-2),1)
            r3=substr($4,(i+addto4-1),1)
            r4=substr($4,i+addto4,1)
            r5="."; if (i<num_char) { r5=substr($4,(i+addto4+1),1) }
            a1="."; if (i>3) { a1=substr($5,(i+addto5-3),1) }
            a2=substr($5,(i+addto5-2),1)
            a3=substr($5,(i+addto5-1),1)
            a4=substr($5,i+addto5,1)
            a5="."; if (i<num_char) { a5=substr($5,(i+addto5+1),1) }
            if ((r1==a1)&&(r2!=a2)&&(r3!=a3)&&(r4!=a4)&&(r5==a5)) { # 3 adjacent SNPs makes a 3bp MNV that has 3 SNPs in it
              printf $1 OFS ($2+i-3+addto4) OFS $3 OFS r2 r3 r4 OFS a2 a3 a4; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
            }
            if ((r1==a1)&&(r2!=a2)&&(r3==a3)&&(r4!=a4)&&(r5==a5)) { # SNP then non-SNP then SNP makes a 3bp MNV that has 2 SNPs in it not 3 SNPs
              printf $1 OFS ($2+i-3+addto4) OFS $3 OFS r2 r3 r4 OFS a2 a3 a4; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
            }
          }
        }
      }
    }

    # look for 3bp MNVs regardless of whether nucleotides match the REF on either side
    # only look in platypus variants having same length REF and ALT, so that we do not get a lot of spurious sub-variants from stuttering
    if ((length($4) >= 3) && (length($4) == length($5))) {
      num_char=length($4)
      for (i=3; i<=num_char; i++) {
        r1=substr($4,(i-2),1)
        r2=substr($4,(i-1),1)
        r3=substr($4,i,1)
        a1=substr($5,(i-2),1)
        a2=substr($5,(i-1),1)
        a3=substr($5,i,1)
        if ((r1!=a1)&&(r2!=a2)&&(r3!=a3)) { # 3 adjacent SNPs makes a 3bp MNV that has 3 SNPs in it
          printf $1 OFS ($2+i-3) OFS $3 OFS r1 r2 r3 OFS a1 a2 a3; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
        }
        if ((r1!=a1)&&(r2==a2)&&(r3!=a3)) { # SNP then non-SNP then SNP makes a 3bp MNV that has 2 SNPs in it not 3 SNPs
          printf $1 OFS ($2+i-3) OFS $3 OFS r1 r2 r3 OFS a1 a2 a3; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
        }
      }
    }

    # look for MNVs that are two SNPs that are a distance of 3 to 10 bp apart
    # we already caught dist=1 in 2bp MNVs and and some of the dist=2 in 3bp MNVs
    # these 2 SNPs may be embedded in a platypus variant that is an indel
    # make sure the nucleotides in-between the 2 SNPs are the same (REF and ALT are both REF) because we can annotate with with gnomad MNV frequency
    # currently we cannot annotate gnomad distance MNVs when the in-between nucleotides are not the same, even though gnomad include that in their frequencies
    # thus we will not output those gnomad distance MNVs that we cannot annotate with gnomad frequencies
    # only look from left to right
    if (length($4) >= 3) {
      for (dist=2; dist<10; dist++) { # this "dist" variable contains gnomad "dist" plus 1
        dist_plus_1 = dist + 1
        if ((dist+1) <= num_char) {
          for (i=dist_plus_1; i<=num_char; i++) {
            snp1_ref=substr($4,(i-dist),1)
            snp2_ref=substr($4,i,1)
            snp1_alt=substr($5,(i-dist),1)
            snp2_alt=substr($5,i,1)
            if ((snp1_ref!=snp1_alt)&&(snp2_ref!=snp2_alt)) {
              middle_is_same=1
              for (j=i-dist+1; j<i; j++) {
                ref_bp=substr($4,j,1)
                alt_bp=substr($5,j,1)
                if (ref_bp!=alt_bp) {
                  middle_is_same=0
                }
              }
              if (middle_is_same==1) {
                mnv_ref=substr($4,(i-dist),(dist+1))
                mnv_alt=substr($5,(i-dist),(dist+1))
                printf $1 OFS ($2+i-dist-1) OFS $3 OFS mnv_ref OFS mnv_alt; for(j=6;j<=NF;++j) printf OFS $j; printf "\n"
              }
            }
          }
        }
      }
    }

  }
}

