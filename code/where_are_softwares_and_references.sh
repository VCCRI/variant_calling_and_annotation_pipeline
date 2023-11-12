#!/bin/bash

export root_of_reference_data=/my/directory/variant_calling_and_annotation_pipeline/reference_data
export root_of_References_and_Databases=/g/data/a32/References_and_Databases
export sw=/g/data/jb96/software
export batch=MY_BATCH # an identifier that will be used to prefix file names that are for aggregation files of all samples in this processing run, rather than the sample level
export R_LIBS_USER=/my/directory/variant_calling_and_annotation_pipeline/victorchang_scripts/R_libraries_for_victorchang_scripts

export genome_version=hg38 # hg38 hg19
export limit_num_of_pbs_submissions=0
export ref_fasta="${root_of_References_and_Databases}"/hg38.noalt.decoy.bwa/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
export ref_fasta_fa="${root_of_References_and_Databases}"/hg38.noalt.decoy.bwa/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa
export ref_fasta_with_alts="${root_of_References_and_Databases}"/Homo_sapiens_assembly38.fasta/Homo_sapiens_assembly38.fasta
export ref_fasta_chroms="${root_of_References_and_Databases}"/hg38.noalt.decoy.bwa.chroms
export humandb="${sw}"/annovar/humandb
export gatk_path="${sw}"/GATK/gatk-4.2.0.0/gatk
export gatk3_path="${sw}"/GATK/GenomeAnalysisTK-3.8-1-0/GenomeAnalysisTK.jar
export bwa="${sw}"/bwa-0.7.17/bwa
export picard_jar="${sw}"/picard/picard-2.18.26/picard.jar
# For GATK HaplotypeCaller:
export gatk_dbsnp="${root_of_References_and_Databases}"/GATK_bundle/hg38/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
# For GATK VariantRecalibrator:
export gatk_db_dir="${root_of_References_and_Databases}"/GATK_bundle/hg38
export mills_resource_vcf="${gatk_db_dir}"/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
export axiomPoly_resource_vcf="${gatk_db_dir}"/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
export dbsnp_resource_vcf="${gatk_db_dir}"/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
export hapmap_resource_vcf="$gatk_db_dir/hapmap_3.3.hg38.vcf.gz"
export omni_resource_vcf="$gatk_db_dir/1000G_omni2.5.hg38.vcf.gz"
export one_thousand_genomes_resource_vcf="$gatk_db_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
# For align bam bqsr
export dbsnp="${root_of_References_and_Databases}"/GATK_bundle/hg38/beta/Homo_sapiens_assembly38.dbsnp138.vcf.gz
export known_indels="${root_of_References_and_Databases}"/GATK_bundle/hg38/beta/Homo_sapiens_assembly38.known_indels.vcf.gz
export gold_std_indels="${root_of_References_and_Databases}"/GATK_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# For annovar annotation
export wgEncodeGencodeBasic=wgEncodeGencodeBasicV33
export EHRF_regions=EHRF_r104_regions
export exac03_ed=exac03_ed1_high
# For vep annotation
export vep_assembly=GRCh38
export vep_loftee="${sw}"/vep_plugins/LOFTEE_GRCh38/loftee
export vep_loftee_human_ancestor_fa="${sw}"/vep_plugins/LOFTEE_GRCh38/human_ancestor/human_ancestor.fa.gz
export vep_loftee_conservation_file="${sw}"/vep_plugins/LOFTEE_GRCh38/conservation_file/phylocsf_gerp.sql
export vep_utrannotator="${sw}"/vep_plugins/UTRannotator/UTRannotator/uORF_starts_ends_GRCh38_PUBLIC.txt
# For gridss structural variant calling
export gridss_jar="${sw}"/gridss_more/2.8.3/gridss-2.8.3-gridss-jar-with-dependencies.jar
export gridss_blacklist="${sw}"/gridss_resources/hg38-blacklist.v2.bed
export regions_of_LowComplexity_SimpleRepeats="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_Repeats_20170921_LowComplexity_SimpleRepeats.bed
# For manta structural variant calling
export manta_configManta="${sw}"/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py
export manta_convertInversion="${sw}"/manta/manta-1.6.0.centos6_x86_64/libexec/convertInversion.py
export manta_samtools=/apps/samtools/1.10/bin/samtools
# For cnvnator copy number variant calling
export cnvnator_executable="${sw}"/CNVnator_gadi/CNVnator/cnvnator
export rootsys_for_cnvnator="${sw}"/root_gadi_from_git/root/root-build
# For structural variant annotation and gene list creation
export SVannotation_software_directory="${sw}"/public_domain_MGRB_structural_variant_annotation/structural_variant_annotation
export ucsc_refseq_cdsexons="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_CDSexons_RefSeq_20211116.txt
export ucsc_refseq_genes="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20211116_plus_TAPVR1.txt
export ucsc_refseq_introns="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_introns_RefSeq_20210505.txt
export ucsc_gencode_cdsexons="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_CDSexons_wgEncodeGencodeCompV38_20211116.txt
export ucsc_gencode_genes="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_genes_wgEncodeGencodeCompV38_20211116.txt
export ucsc_gencode_introns="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_introns_wgEncodeGencodeCompV38_20211116.txt
export refdata_EHRFr104="${sw}"/annovar/humandb/hg38_EHRF_r104_chr.txt
export refdata_gnomadsv="${sw}"/annovar/humandb/hg38_gnomAD-SV_short_ordered_chr.txt
export refdata_DGVgoldStdVar_r20160516="${sw}"/annovar/humandb/hg38_DGVgoldStdVars_r20160516.txt
export refdata_segmental_duplications="${sw}"/annovar/humandb/GRCh38GenomicSuperDup.tab
export refdata_dbscSNV="${sw}"/annovar/humandb/hg38_splice_dbscSNV_1_0p6.txt
export refdata_FANTOM5_TSS="${sw}"/annovar/humandb/hg38_FANTOM5_CAGEr_humanHeart_transcriptStartSites_1column.txt
export refdata_FANTOM5_TSS_5cols="${sw}"/annovar/humandb/hg38_FANTOM5_CAGEr_humanHeart_transcriptStartSites_5columns.txt
export pext_artery_aorta="${sw}"/annovar/humandb/hg38_pext_sv_for_exons_of_Artery_Aorta.bed
export pext_artery_coronary="${sw}"/annovar/humandb/hg38_pext_sv_for_exons_of_Artery_Coronary.bed
export pext_heart_atrial_appendage="${sw}"/annovar/humandb/hg38_pext_sv_for_exons_of_Heart_AtrialAppendage.bed
export pext_heart_left_ventricle="${sw}"/annovar/humandb/hg38_pext_sv_for_exons_of_Heart_LeftVentricle.bed
export pext_genes="${sw}"/annovar/humandb/hg38_pext_sv_for_genes.bed
export hotspot_exons="${sw}"/annovar/humandb/hg38_hotspot_exons_sv.bed
# For splicing/splice sites
export pythonpath_for_spliceai="${sw}"/spliceai_pip_install/spliceai
export spliceogen_sw="${sw}"/spliceogen_2020march_updated_2021march
export fasta_for_spliceogen="${root_of_References_and_Databases}"/hg38.fa_faidx/hg38.fa
export gtf_for_spliceogen="${sw}"/spliceogen_2020july/resources/gencode.v33.basic.annotation.gtf.gz
export ese_for_spliceogen="${sw}"/spliceogen_2020july/resources/ESE.txt
export ess_for_spliceogen="${sw}"/spliceogen_2020july/resources/ESS.txt
# For spliceai tensorflow program
export spliceai_reference="${root_of_References_and_Databases}"/grch38.p12/GRCh38.p12.genome.fa
# For ConanVarvar
export conanvarvar_sv_dir=/g/data/jb96/mg8150/test/ConanVarvar
# For liftover
export liftover_chain="${sw}"/liftover/hg38ToHg19.over.chain.gz
export liftover_to_genome=hg19
# For annotation of genes by bedtools intersect of gene regions, has multiple annotation columns
export gnomad_constraints="${sw}"/annovar/humandb/hg38_gnomadv211_lof_metrics_by_gene.txt
# For gene lists
export exons_chrom_start_end_gene_strand="${root_of_reference_data}"/genes_and_chromosome_positions/UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20211116_plus_TAPVR1.txt
# For platypus, to call multi-nucleotide variants (MNV)
export platypus="${sw}"/Platypus_0.8.1_on_gadi/Platypus_0.8.1/Platypus.py
# For verifybamid
export verifybamid_population_frequencies=/g/data/jb96/software/verifybamid2/VerifyBamID/resource/1000g.phase3.100k.b38.vcf.gz.dat
