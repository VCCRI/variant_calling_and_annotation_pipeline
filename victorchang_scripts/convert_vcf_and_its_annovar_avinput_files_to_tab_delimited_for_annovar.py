#!/usr/bin/python
# python3 convert_vcf_and_its_annovar_avinput_files_to_tab_delimited_for_annovar.py -i input_VCF_file -a input_avinput -o output_tsv_file -p PREFIX_

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2020, Victor Chang Medical Research Institute'
# Modified from a public domain Garvan Institute of Medical Research script by same author.

# perl convert2annovar.pl -format vcf4 humandb/hg38_clinvar_20200224.vcf > humandb/hg38_clinvar_20200224_convert2annovar.avinput
# python3 convert_vcf_and_its_annovar_avinput_files_to_tab_delimited_for_annovar.py -i humandb/hg38_clinvar_20200224.vcf -a humandb/hg38_clinvar_20200224_convert2annovar.avinput -o humandb/hg38_clinvar_20200224.txt -p CLINVAR20200224_
# perl index_annovar.pl humandb/hg38_clinvar_20200224.txt -outfile humandb/hg38_clinvar_20200224.txt.idx --skipsort
#
# table_annovar_ed.pl my_sample.vcf humandb/ -vcfinput -buildver hg38 \
#   -out my_sample_annotated -remove \
#   -protocol clinvar_20200224 \
#   -operation f -nastring . \
#   -arg -time

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
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file having INFO fields annotations, and also read in its annovar avinput file, and output this record with INFO fields as tab-separated columns in annovar data format, with first five columns Chr Start End Ref Alt from avinput and the rest of the columns from the vcf.info fields.')
	parser.add_argument('-i', action="store", dest="in_VCF", required=True, help='Input VCF file')
	parser.add_argument('-a', action="store", dest="in_avinput", required=True, help='Input avinput file, created by running perl convert2annovar.pl -format vcf4 humandb/hg38_my_data.vcf > humandb/hg38_my_data_convert2annovar.avinput')
	parser.add_argument('-o', action="store", dest="out_TSV", required=False, help='Output tab-delimited annovar database text file - one line per VCF record')
	parser.add_argument('-p', action="store", dest="out_column_prefix", required=False, help='If present, prefix this to each output column heading (excluding first five columns)')
	args = parser.parse_args()

	# read each record of input VCF file and process it

	in_VCF = open(args.in_VCF, 'r')
	in_avinput = open(args.in_avinput, 'r')
	out_TSV = open(args.out_TSV, 'w')
	out_column_prefix = ''
	if (args.out_column_prefix is not None):
		out_column_prefix = str(args.out_column_prefix)

	in_header = True
	out_TSV_header = "#Chr\tStart\tEnd\tRef\tAlt"
	sample_FORMAT_part_of_out_TSV_header = ""

	INFO_field_order = []
	INFO_field_count_of_subfields = {}
	INFO_field_names_of_subfields = {}

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
								out_TSV_header = out_TSV_header + "\t" + out_column_prefix + INFO_header
								INFO_field_order.append( INFO_header )
								INFO_field_count_of_subfields[INFO_header] = 1
							else:
								info_subfield_keys = []
								for idx0 in range( 0, num_output_subfields ):
									idx1 = idx0 + 1
									INFO_subfield_name = INFO_header + '_' + str(idx1)
									out_TSV_header = out_TSV_header + "\t" + out_column_prefix + INFO_subfield_name
									INFO_field_order.append( INFO_subfield_name )
									info_subfield_keys.append( INFO_subfield_name )
								INFO_field_count_of_subfields[INFO_header] = num_output_subfields
								INFO_field_names_of_subfields[INFO_header] = info_subfield_keys
						else:

							# Only one subfield column will be output for this INFO field, regardless of how many each variants has - can vary for each allele

							out_TSV_header = out_TSV_header + "\t" + out_column_prefix + INFO_header
							INFO_field_order.append( INFO_header )
							INFO_field_count_of_subfields[INFO_header] = 1

					elif (inline[0:6] == '#CHROM'):

						infields = inline.split("\t")
						num_samples = len(infields) - 9
						if (num_samples < 0):
							num_samples = 0
						if (len(infields) > 9):
							for i in range( 9, len(infields) ):
								out_TSV_header = out_TSV_header + "\t" + out_column_prefix + str(infields[i])

			else:

				# when finished reading in VCF header records, write out the Excel header record

				in_header = False

				out_TSV.write(out_TSV_header + "\n")

		if (in_header == False):

			inline2 = in_avinput.readline()
			infields = inline2.split("\t")
			chrom2 = str(infields[0])
			start2 = str(infields[1])
			end2 = str(infields[2])
			ref2 = str(infields[3])
			alt2 = str(infields[4])

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

			# outline = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs
			outline = chrom2 + "\t" + start2 + "\t" + end2 + "\t" + ref2 + "\t" + alt2

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

			# INFO_field_count_of_subfields[INFO_header]

			for i in range( 0, len(INFO_field_order) ): # For every info field in header, write out this variant's value for it
				this_info_key = INFO_field_order[i]
				this_info_value = '.'
				if this_info_key in INFO_fields_and_subfields:
					this_info_value = str(INFO_fields_and_subfields[this_info_key])
				outline = outline + "\t" + this_info_value

			if (len(infields) > 9) : # For every sample, write out its first field, which is probably the GT field
				for i in range( 9, len(infields) ):
					this_sample = str(infields[i])
					this_sample_fields = this_sample.split(":")
					this_sample_GT = str(this_sample_fields[0])
					outline = outline + "\t" + this_sample_GT

			out_TSV.write(outline + "\n")

	out_TSV.close()


if __name__=='__main__':
    main()

