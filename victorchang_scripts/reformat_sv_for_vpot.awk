# Author: Emma Rath
# Date: 2020-05-20
#
# awk -v sample=19W001486 -f temp.awk 19W001486_gridss_annotated.tsv > 19W001486_gridss_annotated.reformatted_for_vpot.tsv
#
# Reformat an annotated SV or CNV file for input into VPOT.
#
# The output needs to have the following:
# Header line to start with:   #CHROM   POS   END   REF   ALT
# The last 2 columns need to be:   FORMAT   individual_sample_ID
# The FORMAT column of each data line to contain:   GT:DP
# The sample_ID depth component needs to be 10 or more.
# There needs to be a max_gnomadsv_af field containing the highest gnomad AF from gnomadsv_SVTYPE_AF field.
# Thus, this program adds on the following 3 columns:   max_gnomadsv_af   FORMAT   individual_sample_ID

# GRIDSS input file:
# chrom   Start      End        Ref   Alt       Qual     SVEND      SVLEN   SVTYPE  ...   gridssREF   gridssREFPAIR   ...   genotype   ...   gnomadsv_SVTYPE_AF       ...
# chr1    21254075   21254075   A     <INDEL>   1241.9   21254096   0       INDEL   ...   37,40       13,9            ...   0/1        ...   DEL_0.8220340013504028   ...
# chr10   66524187   66524189   T     <DUP:INS> 968.47   66524189   41      DUP     ...   13,27       19,21           ...   0/1        ...   DEL_0.00014000000373926014, DUP_0.05263200029730797   ...

# GRIDSS VAF calculation:
#                ##INFO=<ID=AS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint">
#                ##INFO=<ID=IC,Number=1,Type=Integer,Description="Count of read indels supporting breakpoint">
#                ##INFO=<ID=RAS,Number=1,Type=Integer,Description="Count of assemblies supporting breakpoint from remote breakend">
#                ##INFO=<ID=REF,Number=1,Type=Integer,Description="Count of reads mapping across this breakend">
#                ##INFO=<ID=REFPAIR,Number=1,Type=Integer,Description="Count of reference read pairs spanning this breakpoint supporting the reference allele">
#                ##INFO=<ID=RP,Number=1,Type=Integer,Description="Count of read pairs supporting breakpoint">
#                ##INFO=<ID=SR,Number=1,Type=Integer,Description="Count of split reads supporting breakpoint">
# new_BNDVAF_denominator = this_REFPAIR + this_SR + this_RP + this_IC + this_AS
# new_BNDVAF = float(this_SR + this_RP + this_IC + this_AS) / float(new_BNDVAF_denominator)
# Thus the read depth is gridssREF + gridssREFPAIR. There are 2 values of each, one for each breakend contriuting to the SV. For INV there are 4.
# Add 1st gridssREF to 1st gridssREFPAIR. Add 2nd gridssREF to 2nd gridssREPAIR. Take highest as the depth.

function max(a, b) {
    if (a>b) {
        return a
    } else {
        return b
    }
}

function abs(a) {
    abs_a = a
    if (a < 0) {
        abs_a = a * -1
    }
    return abs_a
}

function extract_depth(SVTYPE, gridssREF, gridssREFPAIR) {
    split(gridssREF,a,",")
    split(gridssREFPAIR,b,",")
    d1=a[1]+b[1]
    d2=a[2]+b[2]
    depth=max(d1,d2)
    if (SVTYPE=="INV") {
        d3=a[3]+b[3]
        d4=a[4]+b[4]
        depth2=max(d3,d4)
        depth=max(depth,depth2)
    }
    if (depth==0) {
        depth=1
    }
    return depth
}

function extract_max_gnomadsv_af(local_gnomadsv_SVTYPE_AF_string) {
    # Two type of input:
    # DEL_0.015557000413537025, DEL_0.09683399647474289, DUP_0.001221999991685152
    # MCNV_0.0009309999877586961,0.160179004073143,0.1831810027360916,0.48351600766181946
    local_highest_af = 0
    if ((local_gnomadsv_SVTYPE_AF_string != ".") && (local_gnomadsv_SVTYPE_AF_string != "")) {
        local_gnomadsv_SVTYPE_AF_string2 = gsub(" ","",local_gnomadsv_SVTYPE_AF_string)
        split(local_gnomadsv_SVTYPE_AF_string, local_gnomadsv_SVTYPE_AF_array, ",")
        for (i=1; i in local_gnomadsv_SVTYPE_AF_array; i++) { 
            local_gnomadsv_SVTYPE_AF = local_gnomadsv_SVTYPE_AF_array[i]
            split(local_gnomadsv_SVTYPE_AF, local_gnomadsv_bits, "_")
            local_gnomadsv_AF = local_gnomadsv_bits[2]
            if (local_gnomadsv_AF == "") {
                local_gnomadsv_AF = local_gnomadsv_bits[1]
            }
            if (local_highest_af == 0) {
                    local_highest_af = local_gnomadsv_AF
                    rtn = local_gnomadsv_AF
            } else {
                if (local_gnomadsv_AF > local_highest_af) {
                    local_highest_af = local_gnomadsv_AF
                }
            }
        }
    }
    return local_highest_af
}

function look_for_Promoter(ensemble_human_regulatory_features_string) {
    hits_Promoter = "NO"
    split(ensemble_human_regulatory_features_string, ensemble_human_regulatory_features_array, ", ")
    for(i=1; i in ensemble_human_regulatory_features_array; i++) {
        if (ensemble_human_regulatory_features_array[i] == "Promoter") {
            hits_Promoter = "YES"
        }
    }
    return hits_Promoter
}

function look_for_hits_in_1_column(col1) {
    hits_columns = "NO"
    if ((col1 == ".") || (col1 == ".")) {
        hits_columns = "NO"
    } else {
        hits_columns = "YES"
    }
    return hits_columns
}

function look_for_hits_in_2_columns(col1, col2) {
    hits_columns = "NO"
    if ((col1 != ".") && (col1 != "")) {
        hits_columns = "YES"
    }
    if ((col2 != ".") && (col2 != "")) {
        hits_columns = "YES"
    }
    return hits_columns
}

BEGIN {
    FS = "\t"
    OFS = "\t"
}
NR==1 {0	YES
    # print the header
    genotype=-1
    SVTYPE=-1
    gridssREF=-1
    gridssREFPAIR=-1
    gnomadsv_SVTYPE_AF=-1
    normalized_read_depth=-1
    ensemble_human_regulatory_features=-1
    segmental_duplication_region_hit_by_left_BND=-1
    segmental_duplication_region_hit_by_right_BND=-1
    exon_containing_left_BND=-1
    exon_containing_right_BND=-1
    gene_CDSexons=-1
    gene_entireGene=-1
    gencode_exon_containing_left_BND=-1
    gencode_exon_containing_right_BND=-1
    gencode_gene_CDSexons=-1
    gencode_gene_entireGene=-1
    dgv_cnvtype_study=-1
    dbscSNV=-1
    otherSamples_exactMatch=-1
    otherSamples_90percentMatch=-1
    includeLowQual_otherSamples_exactMatch=-1
    includeLowQual_otherSamples_90percentMatch=-1
    cnv_size=-1
    SVLEN=-1
    intron_too_short=-1
    gencode_intron_too_short=-1
    for (i=1; i<=NF; ++i) {
        if ($i=="genotype") { genotype=i }
        if ($i=="SVTYPE") { SVTYPE=i }
        if ($i=="gridssREF") { gridssREF=i }
        if ($i=="REF") { gridssREF=i }
        if ($i=="INFO_REF") { gridssREF=i }
        if ($i=="gridssREFPAIR") { gridssREFPAIR=i }
        if ($i=="REFPAIR") { gridssREFPAIR=i }
        if ($i=="INFO_REFPAIR") { gridssREFPAIR=i }
        if ($i=="gnomadsv_SVTYPE_AF") { gnomadsv_SVTYPE_AF=i }
        if ($i=="normalized_read_depth") { normalized_read_depth=i }
        if ($i=="ensemble_human_regulatory_features") { ensemble_human_regulatory_features=i }
        if ($i=="segmental_duplication_region_hit_by_left_BND") { segmental_duplication_region_hit_by_left_BND=i }
        if ($i=="segmental_duplication_region_hit_by_right_BND") { segmental_duplication_region_hit_by_right_BND=i }
        if ($i=="exon_containing_left_BND") { exon_containing_left_BND=i }
        if ($i=="exon_containing_right_BND") { exon_containing_right_BND=i }
        if ($i=="gene_CDSexons") { gene_CDSexons=i }
        if ($i=="gene_entireGene") { gene_entireGene=i }
        if ($i=="gencode_exon_containing_left_BND") { gencode_exon_containing_left_BND=i }
        if ($i=="gencode_exon_containing_right_BND") { gencode_exon_containing_right_BND=i }
        if ($i=="gencode_gene_CDSexons") { gencode_gene_CDSexons=i }
        if ($i=="gencode_gene_entireGene") { gencode_gene_entireGene=i }
        if ($i=="dgv_cnvtype_study") { dgv_cnvtype_study=i }
        if ($i=="dbscSNV") { dbscSNV=i }
        if ($i=="otherSamples_exactMatch") { otherSamples_exactMatch=i }
        if ($i=="otherSamples_90percentMatch") { otherSamples_90percentMatch=i }
        if ($i=="includeLowQual_otherSamples_exactMatch") { includeLowQual_otherSamples_exactMatch=i }
        if ($i=="includeLowQual_otherSamples_90percentMatch") { includeLowQual_otherSamples_90percentMatch=i }
        if ($i=="cnv_size") { cnv_size=i }
        if ($i=="SVLEN") { SVLEN=i }
        if ($i=="svlen") { SVLEN=i }
        if ($i=="CNV_SIZE") { SVLEN=i }
        if ($i=="cnv_size") { SVLEN=i }
        if ($i=="intron_too_short") { intron_too_short=i }
        if ($i=="gencode_intron_too_short") { gencode_intron_too_short=i }
    }
    # printf "#CHROM" FS "POS" FS "END" FS "REF" FS "ALT"
    # for (i=6; i<=NF; ++i) printf FS $i
    printf "#CHROM" FS "POS" FS "END"
    for (i=4; i<=NF; ++i) printf FS $i
    # gridss may or may not have includeLowQual_otherSamples_exactMatch and includeLowQual_otherSamples_90percentMatch. manta and cnvnator wont have it.
    if ((includeLowQual_otherSamples_exactMatch != -1) && (includeLowQual_otherSamples_90percentMatch != -1)) {
      print FS "max_gnomadsv_af" FS "overlaps_promoter" FS "BND_hits_segmental_duplication_region" FS "BND_hits_exon" FS "overlaps_CDSexon" FS "overlaps_gene" FS "overlaps_dgv" FS "overlaps_dbscSNV" FS "has_intron_too_short" FS "overlaps_otherSamples_exactMatch" FS "overlaps_otherSamples_90percentMatch" FS "overlaps_includeLowQual_otherSamples_exactMatch" FS "overlaps_includeLowQual_otherSamples_90percentMatch" FS "any_exactMatch_overlap" FS "any_90percentMatch_overlap" FS "abs_svlen" FS "FORMAT" FS sample
    } else {
      print FS "max_gnomadsv_af" FS "overlaps_promoter" FS "BND_hits_segmental_duplication_region" FS "BND_hits_exon" FS "overlaps_CDSexon" FS "overlaps_gene" FS "overlaps_dgv" FS "overlaps_dbscSNV" FS "has_intron_too_short" FS "overlaps_otherSamples_exactMatch" FS "overlaps_otherSamples_90percentMatch" FS "abs_svlen" FS "FORMAT" FS sample
    }
    next
}
NR>1 {
    this_genotype = "0/1" # Make the default genotype a value, not "./.", because VPOT rejects ./. records
    # get genotype for gridss
    if (genotype != -1) {
        this_genotype = $genotype
    }
    # get genotype for cnvnator
    if (normalized_read_depth != -1) {
        if ($normalized_read_depth <= 0.25) {
            this_genotype = "1/1"
        }
        if ($normalized_read_depth >= 1.25) {
            this_genotype = "0/1"
        }
        if ($normalized_read_depth >= 1.25) {
            this_genotype = "1/1"
        }
    }
    depth=1
    # get depth for gridss
    if ((SVTYPE != -1) && (gridssREF != -1) && (gridssREFPAIR != -1)) {
        depth=extract_depth($SVTYPE, $gridssREF, $gridssREFPAIR)
    }
    if (depth==0) {
        depth = 1
    }
    # fill format fields with depth and genotype
    output_FORMAT="GT:DP"
    output_sample_field=this_genotype":"depth
    # gnomad sv
    output_max_gnomadsv_af = "."
    if (gnomadsv_SVTYPE_AF != -1) {
        output_max_gnomadsv_af = extract_max_gnomadsv_af($gnomadsv_SVTYPE_AF)
    }
    output_overlaps_promoter = "NO"
    if (ensemble_human_regulatory_features != -1) {
        output_overlaps_promoter = look_for_Promoter($ensemble_human_regulatory_features)
    }
    output_BND_hits_segmental_duplication_region = "NO"
    if ((segmental_duplication_region_hit_by_left_BND != -1) && (segmental_duplication_region_hit_by_right_BND != -1)) {
        output_BND_hits_segmental_duplication_region = look_for_hits_in_2_columns($segmental_duplication_region_hit_by_left_BND, $segmental_duplication_region_hit_by_right_BND)
    }
    output_BND_hits_exon = "NO"
    if ((exon_containing_left_BND != -1) && (exon_containing_right_BND != -1)) {
        output_BND_hits_exon = look_for_hits_in_2_columns($exon_containing_left_BND, $exon_containing_right_BND)
    }
    output_overlaps_CDSexon = "NO"
    if (gene_CDSexons != -1) {
        output_overlaps_CDSexon = look_for_hits_in_1_column($gene_CDSexons)
    }
    output_overlaps_gene = "NO"
    if (gene_entireGene != -1) {
        output_overlaps_gene = look_for_hits_in_1_column($gene_entireGene)
    }
    output_BND_hits_gencode_exon = "NO"
    if ((gencode_exon_containing_left_BND != -1) && (gencode_exon_containing_right_BND != -1)) {
        output_BND_hits_gencode_exon = look_for_hits_in_2_columns($gencode_exon_containing_left_BND, $gencode_exon_containing_right_BND)
    }
    output_overlaps_gencode_CDSexon = "NO"
    if (gencode_gene_CDSexons != -1) {
        output_overlaps_gencode_CDSexon = look_for_hits_in_1_column($gencode_gene_CDSexons)
    }
    output_overlaps_gencode_gene = "NO"
    if (gencode_gene_entireGene != -1) {
        output_overlaps_geencode_gene = look_for_hits_in_1_column($gencode_gene_entireGene)
    }
    output_overlaps_dgv = "NO"
    if (dgv_cnvtype_study != -1) {
        output_overlaps_dgv = look_for_hits_in_1_column($dgv_cnvtype_study)
    }
    output_overlaps_dbscSNV = "NO"
    if (dbscSNV != -1) {
        output_overlaps_dbscSNV = look_for_hits_in_1_column($dbscSNV)
    }
    output_has_intron_too_short = "NO"
    if (intron_too_short != -1) {
      if (($intron_too_short != ".") && ($intron_too_short != "")) {
        output_has_intron_too_short = "YES"
      }
    }
    if (gencode_intron_too_short != -1) {
      if (($gencode_intron_too_short != ".") && ($gencode_intron_too_short != "")) {
        output_has_intron_too_short = "YES"
      }
    }
    output_overlaps_otherSamples_exactMatch = "NO"
    if (otherSamples_exactMatch != -1) {
        output_overlaps_otherSamples_exactMatch = look_for_hits_in_1_column($otherSamples_exactMatch)
    }
    output_overlaps_otherSamples_90percentMatch = "NO"
    if (otherSamples_90percentMatch != -1) {
        output_overlaps_otherSamples_90percentMatch = look_for_hits_in_1_column($otherSamples_90percentMatch)
    }
    output_overlaps_includeLowQual_otherSamples_exactMatch = "NO"
    if (includeLowQual_otherSamples_exactMatch != -1) {
        output_overlaps_includeLowQual_otherSamples_exactMatch = look_for_hits_in_1_column($includeLowQual_otherSamples_exactMatch)
    }
    output_overlaps_includeLowQual_otherSamples_90percentMatch = "NO"
    if (includeLowQual_otherSamples_90percentMatch != -1) {
        output_overlaps_includeLowQual_otherSamples_90percentMatch = look_for_hits_in_1_column($includeLowQual_otherSamples_90percentMatch)
    }
    output_abs_svlen = 0
    if (cnv_size != -1) {
        output_abs_svlen = abs($cnv_size)
    }
    if (SVLEN != -1) {
        output_abs_svlen = abs($SVLEN)
    }
    if ((includeLowQual_otherSamples_exactMatch != -1) && (includeLowQual_otherSamples_90percentMatch != -1)) {
      output_any_exactMatch_overlap = "NO"
      if ((output_overlaps_otherSamples_exactMatch == "YES") || (output_overlaps_includeLowQual_otherSamples_exactMatch == "YES")) {
        output_any_exactMatch_overlap = "YES"
      }
      output_any_90percentMatch_overlap = "NO"
      if ((output_overlaps_otherSamples_90percentMatch == "YES") || (output_overlaps_includeLowQual_otherSamples_90percentMatch == "YES")) {
        output_any_90percentMatch_overlap = "YES"
      }
    }
    # output the new fields with the existing line
    if ((includeLowQual_otherSamples_exactMatch != -1) && (includeLowQual_otherSamples_90percentMatch != -1)) {
      print $0, output_max_gnomadsv_af, output_overlaps_promoter, output_BND_hits_segmental_duplication_region, output_BND_hits_exon, output_overlaps_CDSexon, output_overlaps_gene, output_overlaps_dgv, output_overlaps_dbscSNV, output_has_intron_too_short, output_overlaps_otherSamples_exactMatch, output_overlaps_otherSamples_90percentMatch, output_overlaps_includeLowQual_otherSamples_exactMatch, output_overlaps_includeLowQual_otherSamples_90percentMatch, output_any_exactMatch_overlap, output_any_90percentMatch_overlap, output_abs_svlen, output_FORMAT, output_sample_field
    } else {
      print $0, output_max_gnomadsv_af, output_overlaps_promoter, output_BND_hits_segmental_duplication_region, output_BND_hits_exon, output_overlaps_CDSexon, output_overlaps_gene, output_overlaps_dgv, output_overlaps_dbscSNV, output_has_intron_too_short, output_overlaps_otherSamples_exactMatch, output_overlaps_otherSamples_90percentMatch, output_abs_svlen, output_FORMAT, output_sample_field
    }
}


