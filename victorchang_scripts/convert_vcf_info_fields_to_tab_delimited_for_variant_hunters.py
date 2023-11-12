#!/usr/bin/python
# python3 convert_vcf_info_fields_to_tab_delimited_for_variant_hunters.py -i input_VCF_file -o output_tsv_file

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2020, Victor Chang Medical Research Institute'
# Modified from a public domain Garvan Institute of Medical Research script by same author.

import sys
import os
import argparse
import re

DELIMITER_FOR_INFO_MEINFO_FIELD = ','

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def calculate_vaf_from_DP_and_AD( FORMAT_DP, FORMAT_AD ):

	# FORMAT_AD	FORMAT_DP
	# 6,3		9

	VAF = ''
	FORMAT_AD_1 = ''
	FORMAT_AD_2 = ''
	if (is_integer(FORMAT_DP)):
		FORMAT_DP_int = int(FORMAT_DP)
		FORMAT_AD_array = FORMAT_AD.split(',')
		if (len(FORMAT_AD_array) >= 2):
			if (is_integer(FORMAT_AD_array[0]) and is_integer(FORMAT_AD_array[1])):
				FORMAT_AD_1 = FORMAT_AD_array[0]
				FORMAT_AD_2 = FORMAT_AD_array[1]
				ref_allele_depth = int(FORMAT_AD_1)
				alt_allele_depth = int(FORMAT_AD_2)
				VAF = '0'
				ref_plus_alt_allele_depth = ref_allele_depth + alt_allele_depth
				if (ref_plus_alt_allele_depth > 0):
					VAF = str(alt_allele_depth / ref_plus_alt_allele_depth)

	return VAF, FORMAT_AD_1, FORMAT_AD_2

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file having INFO fields annotations, and output this record with INFO fields as tab-separated columns in annovar data format, with first five columns: Chr Start End Ref Alt.')
	parser.add_argument('-i', action="store", dest="in_VCF", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="out_TSV", required=False, help='Output tab-delimited annovar database text file - one line per VCF record')
	parser.add_argument('-p', action="store", dest="out_column_prefix", required=False, help='If present, prefix this to each output column heading (excluding first five columns)')
	parser.add_argument('-id', action="store", dest="output_snp_id", required=False, help='If YES then 3rd output column is vcf SNP_ID, otherwise it is the end position (default).')
	parser.add_argument('-end_id', action="store", dest="output_end_and_snp_id", required=False, help='If YES then 3rd output column is END and 4th output column is vcf SNP_ID, otherwise it is the end position (default).')
	args = parser.parse_args()

	# read each record of input VCF file and process it

	in_VCF = open(args.in_VCF, 'r')
	out_TSV = open(args.out_TSV, 'w')
	out_column_prefix = ''
	if (args.out_column_prefix is not None):
		out_column_prefix = str(args.out_column_prefix)
	output_snp_id = False
	if (args.output_snp_id is not None):
		args_output_snp_id = args.output_snp_id
		args_output_snp_id = args_output_snp_id.upper()
		if (args_output_snp_id == 'YES'):
			output_snp_id = True
	output_end_and_snp_id = False
	if (args.output_end_and_snp_id is not None):
		args_output_end_and_snp_id = args.output_end_and_snp_id
		args_output_end_and_snp_id = args_output_end_and_snp_id.upper()
		if (args_output_end_and_snp_id == 'YES'):
			output_end_and_snp_id = True

	in_header = True
	out_TSV_header = "CHROM\tPOS\tEND\tREF\tALT\tQUAL\tFILTER"
	if (output_snp_id):
		out_TSV_header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER"
	if (output_end_and_snp_id):
		out_TSV_header = "CHROM\tPOS\tEND\tID\tREF\tALT\tQUAL\tFILTER"
	sample_FORMAT_part_of_out_TSV_header = ""

	INFO_field_order = []
	INFO_field_count_of_subfields = {}
	INFO_field_names_of_subfields = {}

	FORMAT_field_order = []
	FORMAT_field_count_of_subfields = {}
	FORMAT_field_names_of_subfields = {}

	can_calculate_vaf = False

	for inline in in_VCF:

		inline = inline.strip()

		if (in_header == True):

			# process each header line read in

			if (inline[0:1] == '#'):
				in_header = True
				if (len(inline) > 11):

					# Process the INFO header lines.

					if (inline[0:11] == '##INFO=<ID='):

						# Add this INFO field and it's subfields to the Excel header lines

						# ##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
						# ##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">
						# ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
						# ##INFO=<ID=TARGETSITEDUPL,Number=1,Type=String,Description="Did a small part of this recognition site get duplicated or deleted during the insertion event. Values:duplication,deletion,noTSD,unknown">
						# ##INFO=<ID=GenomicRegionGene,Number=4,Type=String,Description="Mobile element info of the form TranscriptId|RegionType|GeneName|PfamId. Multiple PfamId separated by & (And sign). Multiple GenomicRegionGene for same ALT allele separated by exclamation mark !.">

						info_regex = re.compile('##INFO=<ID=(.*)>')
						this_INFO_definition = info_regex.search(inline).group(1)

						bits_by_comma = this_INFO_definition.split(',')
						INFO_header = bits_by_comma[0]
						INFO_Number_field = bits_by_comma[1]
						INFO_Number_field = INFO_Number_field.replace('Number=','')
						INFO_field_subfields = []

						if (is_integer(INFO_Number_field)):

							# This INFO field has a fixed number of subfields. There will be multiple output columns, one column per subfield.

							num_output_subfields = int(INFO_Number_field)
							if (num_output_subfields == 0):
								num_output_subfields = 1
							if (num_output_subfields == 1):
								INFO_field_order.append( INFO_header )
								INFO_field_count_of_subfields[INFO_header] = 1
							else:
								info_subfield_keys = []
								for idx0 in range( 0, num_output_subfields ):
									idx1 = idx0 + 1
									INFO_subfield_name = INFO_header + '_' + str(idx1)
									INFO_field_order.append( INFO_subfield_name )
									info_subfield_keys.append( INFO_subfield_name )
								INFO_field_count_of_subfields[INFO_header] = num_output_subfields
								INFO_field_names_of_subfields[INFO_header] = info_subfield_keys
						else:

							# Only one subfield column will be output for this INFO field, regardless of how many each variants has - can vary for each allele

							INFO_field_order.append( INFO_header )
							INFO_field_count_of_subfields[INFO_header] = 1

					# Process the INFO header lines.

					elif (inline[0:11] == '##FORMAT=<I'):

						# Add this FORMAT fieldto the Excel header lines

						# ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
						# ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">
						# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
						##FORMAT=<ID=GP,Number=G,Type=Float,Description="Phred-scaled posterior probabilities for genotypes as defined in the VCF specification">
						##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
						##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

						# Only one subfield column will be output for this FORMAT field, regardless of how many each variants has - can vary for each allele

						format_regex = re.compile('##FORMAT=<ID=(.*)>')
						this_FORMAT_definition = format_regex.search(inline).group(1)

						bits_by_comma = this_FORMAT_definition.split(',')
						FORMAT_header = str(bits_by_comma[0])

						FORMAT_field_order.append( FORMAT_header )
						FORMAT_field_count_of_subfields[FORMAT_header] = 1

					elif (inline[0:6] == '#CHROM'):

						infields = inline.split("\t")
						num_samples = len(infields) - 9
						if (num_samples < 0):
							num_samples = 0

			else:

				# when finished reading in VCF header records, write out the Excel header record

				in_header = False

				if (('AD' in FORMAT_field_order) and ('DP' in FORMAT_field_order)):
					can_calculate_vaf = True

				for i in range( 0, len(INFO_field_order) ): # For every INFO field in header, print out the info header
					out_TSV_header = out_TSV_header + "\t" + out_column_prefix + 'INFO_' + str(INFO_field_order[i])
				for i in range( 0, len(FORMAT_field_order) ): # For every FORMAT field in header, print out the format header
					out_TSV_header = out_TSV_header + "\t" + out_column_prefix + 'FORMAT_' + str(FORMAT_field_order[i])

				if (can_calculate_vaf):
					out_TSV_header = out_TSV_header + "\t" + 'VAF' + "\t" + 'FORMAT_AD_1_ref_allele_depth' + "\t" + 'FORMAT_AD_2_alt_allele_depth'

				out_TSV.write(out_TSV_header + "\n")

		if (in_header == False):

			# Process each data line read in.
			# For each VCF record, process the multiple ALTs separately, one at a time.
			# This means process the line once for the first ALT, then reprocess the line a second time for the second ALT, etc.
			# For the output file that has one line per VCF record and not one line per ALT (ie. multiple ALTs at same POS will produce 1 output line),
			# then output only the first of each INFO annotation (for some INFO fields, there is only 1 INFO annotation, for others there is one per ALT).
			# For INFO fields that can have more than one set of entries, separated by !, output only the first one.

			# 1       10281   .       A       <INS:ME:L1>     16      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-84|,27|,.;TARGETSITEDUPL=duplication      GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.:.:.:.:.:.:.:.:.
			# 1       10450   .       T       <INS:ME:L1>     15      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-176|710|.;TARGETSITEDUPL=duplication    GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.:.:.:.:.:.:.:.:.
			# 1       1243933 .       G       <INS:ME:ALU>,<INS:ME:L1>        18      .       IMPRECISE;SVTYPE=INS;MEINFO=ALU|-16|17|.;TARGETSITEDUPL=duplication;GenomicRegionGene=ENST00000472541|retained_intron|ACAP3|!ENST00000379031|protein_coding|PUSL1|PF01416,ENST00000472541|retained_
			# 1       1244069 .       G       <INS:ME:L1>,<INS:ME:ALU>        13      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-156|157|.;TARGETSITEDUPL=unknown        GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3SQ   ./.:.:.:.:.:.:.:.:.:.:.:.:.     ./.
			# 1       818025  .       A       <INS:ME:L1>     26      .       IMPRECISE;SVTYPE=INS;MEINFO=L1|-20|20|.;TARGETSITEDUPL=unknown;GenomicRegionGene=ENST00000594233|protein_coding|AL645608.2|     GT:MEI5MB:MEI3MB:MEIDU:MEIPLY:MEI5RB:MEI3RB:MEI5SP:MEI3SP:MEI5PR:MEI3PR:MEI5SQ:MEI3
			# 1       16637573        .       A       <INS:ME:ALU>    25      .       IMPRECISE;SVTYPE=INS;MEINFO=ALU|-135|133|.;TARGETSITEDUPL=deletion;GenomicRegionGene=ENST00000375592|protein_coding|FBXO42|PF01344&PF07646&PF00646!ENST00000478089|processed_transcript|FBXO42|PF00646&PF07

			INFO_field_output = [''] * len(INFO_field_order)

			infields = inline.split("\t")
			CHROM = str(infields[0])
			POS = str(infields[1])
			SNV_ID = str(infields[2])
			REF = str(infields[3])
			all_ALTs = str(infields[4])
			QUAL = str(infields[5])
			FILTER = str(infields[6])

			END = str(int(POS) + len(REF) - 1)

			outline = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL + "\t" + FILTER
			if (output_snp_id):
				outline = CHROM + "\t" + POS + "\t" + SNV_ID + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL + "\t" + FILTER
			if (output_end_and_snp_id):
				outline = CHROM + "\t" + POS + "\t" + END + "\t" + SNV_ID + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL + "\t" + FILTER

			INFO = infields[7]
			INFO_fields_and_subfields = {}
			INFO_fields = INFO.split(';')
			for INFO_field in INFO_fields:
				key_and_value = INFO_field.split('=')
				info_key = key_and_value[0]
				if (len(key_and_value) <= 1):
					info_value = info_key # eg. IMPRECISE key has no value. Value assigned will be the key 'IMPRECISE'
					INFO_fields_and_subfields[info_key] = info_value
				else:
					if (info_key in INFO_field_count_of_subfields): # Ignore this variant's info field if it wasn't in the header
						if (INFO_field_count_of_subfields[info_key] == 1): # eg. SVTYPE has one and only one value. eg. for HITGENES we take only the first value which is for first allele
							info_value = key_and_value[1]
							INFO_fields_and_subfields[info_key] = info_value
						elif (INFO_field_count_of_subfields[info_key] == 0): # for eg. IMPRECISE key that has no value, but should be processed above already, not here
							info_value = info_key # eg. IMPRECISE key has no value. Value assigned will be the key 'IMPRECISE'
							INFO_fields_and_subfields[info_key] = info_value
						else: # eg. MEINFO has 4 subfields. eg. IHOMPOS has 2 subfields
							info_values = key_and_value[1]
							info_subfield_keys = INFO_field_names_of_subfields[info_key]
							info_subfield_values = info_values.split(',')
							for i in range( 0, INFO_field_count_of_subfields[info_key] ):
								info_subfield_key = info_subfield_keys[i]
								info_subfield_value = info_subfield_values[i]
								INFO_fields_and_subfields[info_subfield_key] = info_subfield_value

			for i in range( 0, len(INFO_field_order) ): # For every INFO field in header, write out this variant's value for it
				this_info_key = INFO_field_order[i]
				this_info_value = '.'
				if this_info_key in INFO_fields_and_subfields:
					this_info_value = str(INFO_fields_and_subfields[this_info_key])
				outline = outline + "\t" + this_info_value

			FORMAT_fields_and_subfields = {}
			SAMPLE_fields_and_subfields = {}
			if (len(infields) >= 10):
				FORMAT = infields[8]
				SAMPLE = infields[9]
				FORMAT_fields = FORMAT.split(':')
				SAMPLE_fields = SAMPLE.split(':')
				for i in range( 0, len(FORMAT_fields) ):
					this_FORMAT = str(FORMAT_fields[i])
					this_SAMPLE = str(SAMPLE_fields[i])
					FORMAT_fields_and_subfields[this_FORMAT] = this_SAMPLE

			for i in range( 0, len(FORMAT_field_order) ): # For every FORMAT field in header, write out this SAMPLE's variant's value for it
				this_format_key = FORMAT_field_order[i]
				this_format_value = '.'
				if this_format_key in FORMAT_fields_and_subfields:
					this_format_value = str(FORMAT_fields_and_subfields[this_format_key])
				outline = outline + "\t" + this_format_value

			if (can_calculate_vaf):
				VAF, FORMAT_AD_1, FORMAT_AD_2 = ('.', '.', '.')
				if (('DP' in FORMAT_fields_and_subfields) and ('AD' in FORMAT_fields_and_subfields)):
					VAF, FORMAT_AD_1, FORMAT_AD_2 = calculate_vaf_from_DP_and_AD( FORMAT_fields_and_subfields['DP'], FORMAT_fields_and_subfields['AD'] )
				outline = outline + "\t" + VAF + "\t" + FORMAT_AD_1 + "\t" + FORMAT_AD_2

			out_TSV.write(outline + "\n")

	out_TSV.close()


if __name__=='__main__':
    main()

