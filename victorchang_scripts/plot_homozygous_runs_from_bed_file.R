# module load R/4.0.0
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")

# Plot onto chromosomes the regions specified 

system.file(package="ggplot2")


options(width=250)
library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'sampleid', 'hg38', 'infile.txt', 'outfile.txt' ) # for testing
sample = as.character(args[1])
genome_version = as.character(args[2]) # hg19 or hg38
infile = as.character(args[3])
outfile = as.character(args[4])

sample = as.character(sample)

# hg19 chromosome sizes
hg19_chrom_sizes <- structure(list(V1 = c("chrM", "chr1", "chr2", "chr3", "chr4", 
"chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
"chr20", "chr21", "chr22", "chrX", "chrY"), V2 = c(16571L, 249250621L, 
243199373L, 198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 
146364022L, 141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 
107349540L, 102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 
63025520L, 48129895L, 51304566L, 155270560L, 59373566L)), .Names = c("V1", 
"V2"), class = "data.frame", row.names = c(NA, -25L))

# hg38 chromosome sizes
hg38_chrom_sizes <- structure(list(V1 = c("chrM", "chr1", "chr2", "chr3", "chr4", 
"chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
"chr20", "chr21", "chr22", "chrX", "chrY"), V2 = c(16569L, 248956422L, 
242193529L, 198295559L, 190214555L, 181538259L, 170805979L, 159345973L, 
145138636L, 138394717L, 133797422L, 135086622L, 133275309L, 114364328L, 
107043718L, 101991189L, 90338345L, 83257441L, 80373285L, 58617616L, 
64444167L, 46709983L, 50818468L, 156040895L, 57227415L)), .Names = c("V1", 
"V2"), class = "data.frame", row.names = c(NA, -25L))

# hg19 centromere locations
hg19_centromeres <- structure(list(X.bin = c(23L, 20L, 2L, 1L, 14L, 16L, 1L, 14L, 
1L, 1L, 10L, 1L, 15L, 13L, 1L, 1L, 11L, 13L, 1L, 1L, 1L, 12L, 
10L, 10L), chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", 
"chr6", "chr7", "chr8", "chr9", "chrX", "chrY", "chr10", "chr11", 
"chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
"chr19", "chr20", "chr21", "chr22"), chromStart = c(121535434L, 
92326171L, 90504854L, 49660117L, 46405641L, 58830166L, 58054331L, 
43838887L, 47367679L, 58632012L, 10104553L, 39254935L, 51644205L, 
34856694L, 16000000L, 16000000L, 17000000L, 35335801L, 22263006L, 
15460898L, 24681782L, 26369569L, 11288129L, 13000000L), chromEnd = c(124535434L, 
95326171L, 93504854L, 52660117L, 49405641L, 61830166L, 61054331L, 
46838887L, 50367679L, 61632012L, 13104553L, 42254935L, 54644205L, 
37856694L, 19000000L, 19000000L, 20000000L, 38335801L, 25263006L, 
18460898L, 27681782L, 29369569L, 14288129L, 16000000L), ix = c(1270L, 
770L, 784L, 447L, 452L, 628L, 564L, 376L, 411L, 583L, 105L, 341L, 
447L, 304L, 3L, 3L, 3L, 354L, 192L, 125L, 410L, 275L, 22L, 3L
), n = c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
"N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N"
), size = c(3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L), type = c("centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere"), bridge = c("no", "no", 
"no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", 
"no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no"
)), .Names = c("X.bin", "chrom", "chromStart", "chromEnd", "ix", 
"n", "size", "type", "bridge"), class = "data.frame", row.names = c(NA, 
-24L))

# hg38 centromere locations
hg38_centromeres <- structure(list(X.bin = c(23L, 20L, 2L, 1L, 14L, 16L, 1L, 14L, 
1L, 1L, 10L, 1L, 15L, 13L, 1L, 1L, 11L, 13L, 1L, 1L, 1L, 12L, 
10L, 10L), chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", 
"chr6", "chr7", "chr8", "chr9", "chrX", "chrY", "chr10", "chr11", 
"chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
"chr19", "chr20", "chr21", "chr22"), chromStart = c(122026459L, 
92188145L, 90772458L, 49712061L, 46485900L, 58553888L, 58169653L, 
44033744L, 43389635L, 58605579L, 10316944L, 39686682L, 51078348L, 
34769407L, 16000000L, 16000000L, 17083673L, 36311158L, 22813679L, 
15460899L, 24498980L, 26436232L, 10864560L, 12954788L), chromEnd = c(124932724L, 
94090557L, 93655574L, 51743951L, 50059807L, 59829934L, 61528020L, 
45877265L, 45518558L, 62412542L, 10544039L, 41593521L, 54425074L, 
37185252L, 18051248L, 18173523L, 19725254L, 38265669L, 26616164L, 
20861206L, 27190874L, 30038348L, 12915808L, 15054318L), ix = c(1270L, 
770L, 784L, 447L, 452L, 628L, 564L, 376L, 411L, 583L, 105L, 341L, 
447L, 304L, 3L, 3L, 3L, 354L, 192L, 125L, 410L, 275L, 22L, 3L
), n = c("N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", 
"N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N", "N"
), size = c(3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 3000000L, 
3000000L, 3000000L, 3000000L, 3000000L, 3000000L), type = c("centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere", "centromere", "centromere", 
"centromere", "centromere", "centromere"), bridge = c("no", "no", 
"no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no", 
"no", "no", "no", "no", "no", "no", "no", "no", "no", "no", "no"
)), .Names = c("X.bin", "chrom", "chromStart", "chromEnd", "ix", 
"n", "size", "type", "bridge"), class = "data.frame", row.names = c(NA, 
-24L))

chrom_sizes = hg38_chrom_sizes
centromeres = hg38_centromeres
if (genome_version == "hg19") {
  chrom_sizes = hg19_chrom_sizes
  centromeres = hg19_centromeres
}

# set the column names for the datasets
# IMPORTANT: fields common across datasets should have the same name in each
colnames(chrom_sizes) <- c("chromosome", "size")
colnames(centromeres) <- c('bin', "chromosome", 'start', 'end',
                       'ix', 'n', 'size', 'type', 'bridge')

# create an ordered factor level to use for the chromosomes in all the datasets
chrom_order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                 "chr22", "chrX", "chrY", "chrM")
chrom_key <- setNames(object = as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
                                              12, 13, 14, 15, 16, 17, 18, 19, 20, 
                                              21, 22, 23, 24, 25)), 
                      nm = chrom_order)
chrom_order <- factor(x = chrom_order, levels = rev(chrom_order))


# convert the chromosome column in each dataset to the ordered factor
chrom_sizes[["chromosome"]] <- factor(x = chrom_sizes[["chromosome"]], 
                                      levels = chrom_order)
centromeres[["chromosome"]] <- factor(x = centromeres[["chromosome"]], 
                                      levels = chrom_order)

# Read in the sample data of homozygous run regions as a bed file

data = read.table( infile, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )

colnames(data)[1] = 'homozygosity_run_chrom'
colnames(data)[2] = 'hrun_start'
colnames(data)[3] = 'hrun_end'

if (nrow(data) > 0) {
data$hrun_start = as.numeric(as.character(data$hrun_start))
data$hrun_end = as.numeric(as.character(data$hrun_end))
data$gene = ''
data$feature = 'homozygous run'
names(data)[names(data) == "homozygosity_run_chrom"] = "chromosome"
names(data)[names(data) == "hrun_start"] = "start"
names(data)[names(data) == "hrun_end"] = "end"
data[["chromosome"]] = factor(x = data[["chromosome"]], levels = chrom_order)

ggplot(data = chrom_sizes) + 
    # base rectangles for the chroms, with numeric value for each chrom on the x-axis
    geom_rect(aes(xmin = as.numeric(chromosome) - 0.2, 
                  xmax = as.numeric(chromosome) + 0.2, 
                  ymax = size, ymin = 0), 
              colour="black", fill = "white") + 
    # rotate the plot 90 degrees
    coord_flip() +
    # black & white color theme 
    theme(axis.text.x = element_text(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) + 
    # give the appearance of a discrete axis with chrom labels
    scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
    # add bands for centromeres
    geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2, 
                                      xmax = as.numeric(chromosome) + 0.2, 
                                      ymax = end, ymin = start)) +
    # add bands for feature value
    geom_rect(data = data, aes(xmin = as.numeric(chromosome) - 0.2, 
                                     xmax = as.numeric(chromosome) + 0.2, 
                                     ymax = end, ymin = start, fill = feature)) + 
    ggtitle(sample) +
    # supress scientific notation on the y-axis
    scale_y_continuous(labels = comma) +
    ylab("region (bp)")

ggsave(outfile)

} else {

ggplot(data = chrom_sizes) +
    # base rectangles for the chroms, with numeric value for each chrom on the x-axis
    geom_rect(aes(xmin = as.numeric(chromosome) - 0.2,
                  xmax = as.numeric(chromosome) + 0.2,
                  ymax = size, ymin = 0),
              colour="black", fill = "white") +
    # rotate the plot 90 degrees
    coord_flip() +
    # black & white color theme 
    theme(axis.text.x = element_text(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    # give the appearance of a discrete axis with chrom labels
    scale_x_discrete(name = "chromosome", limits = names(chrom_key)) +
    # add bands for centromeres
    geom_rect(data = centromeres, aes(xmin = as.numeric(chromosome) - 0.2,
                                      xmax = as.numeric(chromosome) + 0.2,
                                      ymax = end, ymin = start)) +
    ggtitle(sample) +
    # supress scientific notation on the y-axis
    scale_y_continuous(labels = comma) +
    ylab("region (bp)")

ggsave(outfile)

}
