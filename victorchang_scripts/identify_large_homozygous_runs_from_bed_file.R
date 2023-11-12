# module load R/4.0.0
# module load intel-compiler/2019.3.199 # may be needed for install.packages("blahblah")

#system.file(package="ggplot2")

options(width=250)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( 'infile.txt', 'outfile.txt' ) # for testing
infile = as.character(args[1])
outfile = as.character(args[2])

data = read.table( infile, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(data)[1] = 'homozygosity_run_chrom'
colnames(data)[2] = 'hrun_start'
colnames(data)[3] = 'hrun_end'

data$homozygosity_run_chrom = as.character(data$homozygosity_run_chrom)
data$hrun_start = as.numeric(as.character(data$hrun_start))
data$hrun_end = as.numeric(as.character(data$hrun_end))
data$hrun_length = data$hrun_end - data$hrun_start

# for each homozygosity run, see whether it is the start of 4,000,000 bp or more of a region having 85% or more homozygosity

# calibrated for AGHA_9, but was too restrictive for CHD_BATCH_9
# min_percent_to_qualify_as_hrun = 85
# min_length_to_qualify_as_hrun = 4000000

# calibrated for CHD_BATCH_9
# min_percent_to_qualify_as_hrun = 70
# min_length_to_qualify_as_hrun = 3000000

#min_percent_to_qualify_as_hrun = 70
#min_length_to_qualify_as_hrun = 3000000

#min_percent_to_qualify_as_hrun = 55
#min_length_to_qualify_as_hrun = 2000000

min_percent_to_qualify_as_hrun = 40
min_length_to_qualify_as_hrun = 1500000

large_hruns = NULL
chrom_list = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
#chrom_list = c("chr14") # debug

for (this_chrom in chrom_list) {
  sset = data[(data$homozygosity_run_chrom==this_chrom),]
##edd  print(c("chromo working -", this_chrom)) # debug
  keep_trying_new_start_positions = 1
  idx_for_this_start = 1
  if (nrow(sset) < 1) {
    keep_trying_new_start_positions = 0
  }

  while (keep_trying_new_start_positions == 1) {

    this_start = sset$hrun_start[idx_for_this_start]
    keep_adding_homozygous_lengths_to_this_start_position = 1
    idx_for_this_end = idx_for_this_start
    prev_end = sset$hrun_end[idx_for_this_start]
    prev_sum_of_homozygous_regions = 0
    prev_percentage_homozygous_in_total_region = 0
    #if (this_start < 100000000) { # debug
    #  keep_adding_homozygous_lengths_to_this_start_position = 0 # debug
    #} # debug
    #if (this_start > 100000000) { # debug
    #  keep_adding_homozygous_lengths_to_this_start_position = 0 # debug
    #  keep_trying_new_start_positions = 0 # debug
    #} # debug

    while (keep_adding_homozygous_lengths_to_this_start_position == 1) {

      this_end = sset$hrun_end[idx_for_this_end]
      this_length = sset$hrun_length[idx_for_this_end]
      sum_of_homozygous_regions = prev_sum_of_homozygous_regions + this_length
      total_length_of_regions = this_end - this_start
      percentage_homozygous_in_total_region = sum_of_homozygous_regions / total_length_of_regions * 100
      #print(c(sset$hrun_start[idx_for_this_end], sset$hrun_end[idx_for_this_end], percentage_homozygous_in_total_region, sum_of_homozygous_regions)) # debug
      if (percentage_homozygous_in_total_region >= min_percent_to_qualify_as_hrun) {
        prev_end = this_end
        prev_sum_of_homozygous_regions = sum_of_homozygous_regions
        prev_percentage_homozygous_in_total_region = percentage_homozygous_in_total_region
      } else {
        keep_adding_homozygous_lengths_to_this_start_position = 0
      }

      #if (keep_adding_homozygous_lengths_to_this_start_position == 1) { # debug
      #  print(c(this_start, "[", sset$hrun_start[idx_for_this_end], sset$hrun_end[idx_for_this_end], data$hrun_length[idx_for_this_end], "]", sum_of_homozygous_regions, total_length_of_regions, percentage_homozygous_in_total_region, "----------------")) # debug
      #} # debug
      idx_for_this_end = idx_for_this_end + 1
      if (idx_for_this_end > nrow(sset)) {
        keep_adding_homozygous_lengths_to_this_start_position = 0
      }
    }

    # We have been accumulating a large homozygous region from smaller ones
    # and now we have hit a smaller run that is too far away from the large run
    # so it is time to write out the large run.
    # But first, adjust the end position.
    if ((keep_adding_homozygous_lengths_to_this_start_position == 0) & (prev_percentage_homozygous_in_total_region >= min_percent_to_qualify_as_hrun) & (prev_sum_of_homozygous_regions >= min_length_to_qualify_as_hrun)) {
      # If we have a very large run, we can overshoot when a small run that really should be considered to be too far away still satisfies the percentage.
      # So let's work backwards and refine the end position of the large homozygous run.
      # If the end smaller runs would not have satisfied the percentage when we accumulate run working backwards, then don't include them.
      true_end = prev_end
      this_end_point_is_good = 0
      keep_looking_for_a_good_end_point = 1
      idx_for_adjusted_end = idx_for_this_end - 1
#
      ##edd print(c("idx_for_adjusted_end", idx_for_adjusted_end, "idx_for_this_end", idx_for_this_end, "this_start", this_start, "idx_for_this_start", idx_for_this_start)) # debug
      if (idx_for_adjusted_end == 1) {     # stop now as idx_for_backwards will be 0 	#edd
         keep_looking_for_a_good_end_point = 0 											#edd
      }                                                                        			#edd
#
      while (keep_looking_for_a_good_end_point == 1) {

        # Let's evaluate this end_point. Going backwards, does it always satisfy the percent_homozygosity criteria?
        #print(c("evaluate true_end", true_end)) # debug
        true_end = sset$hrun_end[idx_for_adjusted_end]
        backwards_sum_of_homozygous_regions = 0
        idx_for_backwards = idx_for_adjusted_end - 1
##edd        print(c("idx_for_backwards", idx_for_backwards, "idx_for_adjusted_end", idx_for_adjusted_end, "idx_for_this_end", idx_for_this_end, "this_start", this_start, "idx_for_this_start", idx_for_this_start)) # debug
        if (idx_for_backwards == this_start) {
          this_end_point_is_good = 1
        } else {
          keep_evaluating_this_end_point = 1
          this_end_point_is_good = 1
          while (keep_evaluating_this_end_point == 1) {
             ##edd print(c("backwards_sum_of_homozygous_regions", backwards_sum_of_homozygous_regions, " sset$hrun_end[idx_for_backwards] - ", sset$hrun_end[idx_for_backwards], "sset$hrun_start[idx_for_backwards] - ", sset$hrun_start[idx_for_backwards])) # debug
             ##edd print(c("true_end", true_end)) # debug
             backwards_sum_of_homozygous_regions = backwards_sum_of_homozygous_regions + (sset$hrun_end[idx_for_backwards] - sset$hrun_start[idx_for_backwards] + 1)
             backwards_total_region = true_end - sset$hrun_start[idx_for_backwards] + 1
             ##edd print(c("backwards_sum_of_homozygous_regions", backwards_sum_of_homozygous_regions)) # debug
             ##edd print(c("backwards_total_region", backwards_total_region)) # debug
             backwards_percent = backwards_sum_of_homozygous_regions / backwards_total_region * 100
             ##edd print(c("problem line - true_end", true_end, "percent", backwards_percent, "backstart", sset$hrun_start[idx_for_backwards], "start", this_start)) # debug
             ##edd print(c("problem line - min_percent_to_qualify_as_hrun", min_percent_to_qualify_as_hrun)) # debug
             if (backwards_percent < min_percent_to_qualify_as_hrun) {
               ##edd print(c("notok:true_end", true_end, "percent", backwards_percent, "backstart", sset$hrun_start[idx_for_backwards], "start", this_start)) # debug
               keep_evaluating_this_end_point = 0
               this_end_point_is_good = 0
             }
             idx_for_backwards = idx_for_backwards - 1
#             print(c("check 2 - idx_for_backwards", idx_for_backwards, "idx_for_this_start", idx_for_this_start)) # debug
             if (idx_for_backwards < idx_for_this_start) {
#               print(c("pass check")) # debug
               keep_evaluating_this_end_point = 0  # stop 
             }
           }
        }

        # We have evaluated this end_point going backwards.
        # If all its percent_homozygosities are good, then we will use it as the end_point for the run_of_homozygosity starting from this_start.
        # If not, then evaluate the end_point prior to it.
        if (this_end_point_is_good == 1) {
          keep_looking_for_a_good_end_point = 0
        } else {
          idx_for_adjusted_end = idx_for_adjusted_end - 1
##edd          print(c("check 3 - idx_for_backwards", idx_for_backwards, "idx_for_adjusted_end", idx_for_adjusted_end)) # debug
#edd          if (idx_for_adjusted_end < idx_for_this_start) {
          if (idx_for_adjusted_end < idx_for_this_start || idx_for_adjusted_end == 1) {     # edd
            keep_looking_for_a_good_end_point = 0
          }
##edd          print(c("check 4 - keep_looking_for_a_good_end_point", keep_looking_for_a_good_end_point)) # debug
        }
      }

      # For this_start and the adjusted end_point, record this run of homozygosity
      if (this_end_point_is_good == 1) {
        #print(c("found a good end point", this_chrom, this_start, true_end)) # debug
        this_large_hrun = as.matrix(t(c( this_chrom, this_start, true_end )))
        large_hruns = rbind( large_hruns, this_large_hrun )
      }
    }

    idx_for_this_start = idx_for_this_end + 1
    if (idx_for_this_start > nrow(sset)) {
      keep_trying_new_start_positions = 0
    }
  }
}
if (is.null(large_hruns)) {

  output_is_empty = "chrom\tstart\tend"
  writeLines( output_is_empty, outfile )

} else {
  large_hruns = as.data.frame(as.matrix(large_hruns))
  colnames(large_hruns) = c("chrom", "start", "end")
  large_hruns$start = as.numeric(as.character(large_hruns$start))
  large_hruns$end = as.numeric(as.character(large_hruns$end))
  write.table( large_hruns, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )
}



