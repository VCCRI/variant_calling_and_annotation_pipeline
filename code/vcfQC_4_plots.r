library(RColorBrewer)
library(tidyr)

setwd('/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/MY_BATCH.plink/')
ped <- '/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/MY_BATCH_ped_file.ped'
inputfile <- '/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/MY_BATCH.plink/MY_BATCH_vcf_istats.txt'
covfile <- '/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/MY_BATCH.vcfQC_3_get_coverage.output.txt'
colours_file <- '/my/directory/variant_calling_and_annotation_pipeline/code/vcfQC_4_plots_MY_BATCH_cols.txt'
pedfile <- read.table(file=ped, header=T)
istats <- read.table(file=inputfile, header=T)
coverage <- read.table(file=covfile, header=T)
family_cols <- read.table(file=colours_file, header=T)
# we have 10 family groups - so only limited palettes with that many colours
cols <- brewer.pal(max(family_cols$colour), 'Set1')
cols_s <- cols[family_cols$colour]
custom_blue <- "#2b8cbe"
custom_red <- "#d7191c"

# Convert sampleID-familyID to just sampleID in case we have a mix of those identifiers in our data
# This happens when VCGS have hyphenated IDs but the pipeline started creating files with non-hyphenated names.
#
pedfile = separate(data=pedfile, col=Subject, into=c("Subject_Sample", "Subject_Family"), sep="-")
names(pedfile)[names(pedfile) == 'Subject_Sample'] <- 'Subject'
pedfile$Subject_Family = NULL
pedfile = separate(data=pedfile, col=Father, into=c("Father_Sample", "Father_Family"), sep="-")
names(pedfile)[names(pedfile) == 'Father_Sample'] <- 'Father'
pedfile$Father_Family = NULL
pedfile = separate(data=pedfile, col=Mother, into=c("Mother_Sample", "Mother_Family"), sep="-")
names(pedfile)[names(pedfile) == 'Mother_Sample'] <- 'Mother'
pedfile$Mother_Family = NULL
#
istats = separate(data=istats, col=ID, into=c("ID_Sample", "ID_Family"), sep="-")
names(istats)[names(istats) == 'ID_Sample'] <- 'ID'
istats$ID_Family = NULL
#
coverage = separate(data=coverage, col=ID, into=c("ID_Sample", "ID_Family"), sep="-")
names(coverage)[names(coverage) == 'ID_Sample'] <- 'ID'
coverage$ID_Family = NULL

#
# Drop any samples = 0 from the pedfile. They were needed for other programs for incomplete trios.
# In this program sample = 0 causes problems with plotting (sample decalage) and factors (pedfile$Subject not same as istats$ID).
new_pedfile = pedfile[(pedfile$Subject != 0),]
pedfile = new_pedfile
pedfile$Subject = factor(pedfile$Subject)

#
stats <- NULL
for (i in 1:nrow(pedfile))
{
  for (j in 1:nrow(istats))
  {
    if (pedfile$Subject[i]==istats$ID[j]) {
      tmp <-cbind(pedfile[i,], istats[j,])
      stats <- rbind(stats, tmp)
      break
    }
  }
}

stats = separate(data=stats, col=Subject, into=c("Subject_Sample", "Subject_Family"), sep="-")
names(stats)[names(stats) == 'Subject_Sample'] <- 'Subject'
stats$Subject_Family = NULL
stats = separate(data=stats, col=Father, into=c("Father_Sample", "Father_Family"), sep="-")
names(stats)[names(stats) == 'Father_Sample'] <- 'Father'
stats$Father_Family = NULL
stats = separate(data=stats, col=Mother, into=c("Mother_Sample", "Mother_Family"), sep="-")
names(stats)[names(stats) == 'Mother_Sample'] <- 'Mother'
stats$Mother_Family = NULL

#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Total_number_of_Variants.png",980,680)
barplot(stats$NALT, cex.name=0.8, cex.axis=0.8, names=stats$ID, las=2, 
        ylab="Number of variants", col=cols_s)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Total_number_of_het_Variants.png",980,680)
barplot(stats$NHET, cex.name=0.8, cex.axis=0.8, names=stats$ID, las=2, 
        ylab="Number of het variants", col=cols_s)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Total_number_of_PASS_Variants.png",980,680)
barplot(stats$PASS, cex.name=0.8, cex.axis=0.8, names=stats$ID, las=2, 
        ylab="Number of variants", col=cols_s)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/TITV_of_Variants.png",980,680)
barplot(stats$TITV, ylim= c(0,2), cex.name=0.8, cex.axis=0.8, 
        names=stats$ID, las=2, ylab="TITV ratio", col=cols_s)
dev.off()
#
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Total_Nbr_compare_PASS_Variants.png",980,680)
barplot(t(stats[c('NALT','PASS')]), beside=T, cex.name=0.8, cex.axis=0.8, 
        names=stats$ID, las=2, ylab="Number of variants", 
        legend.text=c('NALT','PASS'),col=c(custom_blue, custom_red))
dev.off()
#
# coverage plots
#
covstats <- NULL
for (i in 1:nrow(stats))
{
  for (j in 1:nrow(coverage))
  {
    if (stats$Subject[i]==coverage$ID[j]) { #
      tmp <-cbind(stats[i,], coverage[j,])
      covstats <- rbind(covstats, tmp)
      break
    }
  }
}
#
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Tot_vars_vs_avgcov_of_samples.png",980,680)
plot(covstats$NALT, covstats$MEAN_COVERAGE, 
     main="Tot_vars_vs_avgcov_of_samples", xlab="Number of variants", 
     ylab="Coverage", col=cols_s, pch=19, cex=1.5)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/PASS_vars_vs_avgcov_of_samples.png",980,680)
plot(covstats$PASS, covstats$MEAN_COVERAGE, 
     main="PASS_vars_vs_avgcov_of_samples", xlab="Number of PASS variants", 
     ylab="Coverage", col=cols_s, pch=19, cex=1.5)
dev.off()
#
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/HETS_vars_vs_avgcov_of_samples.png",980,680)
plot(covstats$NHET, covstats$MEAN_COVERAGE, 
     main="HET_vars_vs_avgcov_of_samples", xlab="Number of HETS variants", 
     ylab="Coverage", col=cols_s, pch=19, cex=1.5)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Quality_vs_Num_of_Variants_of_samples.png",980,680)
plot(covstats$QUAL, covstats$NALT, main="Quality_vs_Num_of_VARS_for_samples", 
     xlab="Mean QUAL for variants", ylab="Number of variants", col=cols_s, 
     pch=19, cex=1.5)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Quality_vs_Num_of_PASS_Variants_of_samples.png",980,680)
plot(covstats$QUAL, covstats$PASS, 
     main="Quality_vs_Num_of_PASS_vars_for_samples", 
     xlab="Mean QUAL for variants", ylab="Number of PASS variants", col=cols_s, 
     pch=19, cex=1.5)
dev.off()

#
# missingness plots
#
#
missing_file <- "/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/MY_BATCH.plink/MY_BATCH._chrall_missing.imiss"
missing<-read.table(file=missing_file,header=T)
# want output in same order as ped file (consitency with above plots)
# subject id is "IID" in missingness data
n_miss <- missing[pedfile$Subject,]$N_MISS
f_miss <- missing[pedfile$Subject,]$F_MISS
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Nbr_of_missing_genotypes_in_samples.png",1080,680)
barplot(n_miss, cex.name=0.8, cex.axis=0.75, las=2, 
        ylab="Number of missing genotypes", names=pedfile$Subject,
        border="lightblue", col=cols_s)
dev.off()
#
png("/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_qc/Frequency_of_missing_genotypes_in_samples.png",1080,680)
barplot(f_miss, cex.name=0.8, cex.axis=0.75, las=2, 
        ylab="Freq of missing genotypes", names=pedfile$Subject,
        border="lightblue", col=cols_s)
dev.off()
#

