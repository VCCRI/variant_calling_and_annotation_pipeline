#!/usr/bin/python
# python3 convert_vcf_info_fields_to_tab_delimited_for_annovar_qual.py -i input_VCF_file -o output_tsv_file [-add_svlen NO]

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Victor Chang Medical Research Institute'
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
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file having INFO fields annotations, and output this record with INFO fields as tab-separated columns in annovar data format, with first five columns: Chr Start End Ref Alt.')
	parser.add_argument('-i', action="store", dest="in_VCF", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="out_TSV", required=False, help='Output tab-delimited annovar database text file - one line per VCF record')
	parser.add_argument('-p', action="store", dest="out_column_prefix", required=False, help='If present, prefix this to each output column heading (excluding first five columns)')
	parser.add_argument('-add_svlen', action="store", dest="add_svlen", required=False, help='If YES or not present, this program adds SVLEN as the last column')
	args = parser.parse_args()

	# read each record of input VCF file and process it

	in_VCF = open(args.in_VCF, 'r')
	out_TSV = open(args.out_TSV, 'w')
	out_column_prefix = ''
	if (args.out_column_prefix is not None):
		out_column_prefix = str(args.out_column_prefix)
	add_svlen = 'YES'
	if (args.add_svlen is not None):
		add_svlen = args.add_svlen

	in_header = True
	out_TSV_header = "#Chr\tStart\tEnd\tRef\tAlt\tQual"
	sample_FORMAT_part_of_out_TSV_header = ""

	INFO_field_order = []
	INFO_field_count_of_subfields = {}
	INFO_field_names_of_subfields = {}

	for inline in in_VCF:

		#inline = inline.strip()

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
						infields[-1] = infields[-1].strip()
						num_samples = len(infields) - 9
						if (num_samples < 0):
							num_samples = 0
						if (len(infields) > 9): # For every sample, write out the header, which is the sample-id, for its first field, which is probably the GT field
							for i in range( 9, len(infields) ):
								out_TSV_header = out_TSV_header + "\t" + out_column_prefix + str(infields[i])
						out_TSV_header = out_TSV_header + "\t" + out_column_prefix + "VAF_in_INFO" # write out the header for the INFO VAF field, in case there is one
						if (len(infields) > 9): # For every sample, write out the header for its sample_format VAF field, in case there is one
							for i in range( 9, len(infields) ):
								out_TSV_header = out_TSV_header + "\t" + out_column_prefix + "VAF_" + str(infields[i])

			else:

				# When finished reading in VCF header records, write out the Excel header record.

				in_header = False

				svlen_col_name = 'SVLEN'
				if (svlen_col_name in INFO_field_count_of_subfields):
					do_nothing = 1
				else:
					if (add_svlen == 'YES'):
						INFO_field_order.append( svlen_col_name )
						INFO_field_count_of_subfields[ svlen_col_name ] = 1
						out_TSV_header = out_TSV_header + "\t" + out_column_prefix + svlen_col_name

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
			infields[-1] = infields[-1].strip()
			CHROM = str(infields[0])
			POS = str(infields[1])
			SNV_ID = str(infields[2])
			REF = str(infields[3])
			all_ALTs = str(infields[4])
			QUAL = str(infields[5])
			FILTER = str(infields[6])

			END = str(int(POS) + len(REF) - 1)

			outline_part1 = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL
			outline_part2 = ''

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

			# Make any needed adjustments to this VCF record.
			# The following adjustments are needed in manta vcf files:
			#   * fill in SVLEN of INS records that don't have it and do have LEFT_SVINSSEQ and/or RIGHT_SVINSSEQ sequences
			#   * for DUP records that have the same value for POS and END, set the END as POS + SVLEN
			#   * for INV records that have the same value for POS and END, set the END as POS + SVLEN
			#   * for <INS:ME: records, set END to where the TSD (tandem site duplication) ends

			if ('SVTYPE' in INFO_fields_and_subfields):

				if (INFO_fields_and_subfields['SVTYPE'] == 'INS'):
					if ('SVLEN' in INFO_fields_and_subfields):
						do_nothing = 1
					else:
						# Some manta INS records have same value for POS and INFO.END, and no INFO.SVLEN, and have INFO.LEFT_SVINSSEQ and INFO.RIGHT_SVINSSEQ sequences.
						# Get the SVLEN from INFO.LEFT_SVINSSEQ and INFO.RIGHT_SVINSSEQ.
						svlen1 = 0
						svlen2 = 0
						if ('LEFT_SVINSSEQ' in INFO_fields_and_subfields):
							svlen1 = int(len(INFO_fields_and_subfields['LEFT_SVINSSEQ']))
						if ('RIGHT_SVINSSEQ' in INFO_fields_and_subfields):
							svlen2 = int(len(INFO_fields_and_subfields['RIGHT_SVINSSEQ']))
						svlen = max( svlen1, svlen2 )
						if (add_svlen == 'YES'):
							INFO_fields_and_subfields['SVLEN'] = str(svlen)

				if ((INFO_fields_and_subfields['SVTYPE'] == 'DUP') or (INFO_fields_and_subfields['SVTYPE'] == 'INV') or (INFO_fields_and_subfields['SVTYPE'] == 'DEL')):
					if (int(END) <= (int(POS) + 1)): # if POS and END are the same or differ by 1, then we need to find the true end
						if ('END' in INFO_fields_and_subfields):
							END = str(int(INFO_fields_and_subfields['END']))
					if (int(END) <= (int(POS) + 1)): # if END is still the same as POS or is still only 1 more bp, then still need to find the true end
						if ('SVLEN' in INFO_fields_and_subfields):
							this_SVLEN = abs(int(INFO_fields_and_subfields['SVLEN']))
							END = str(int(POS) + this_SVLEN - 1)
					outline_part1 = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL

				if (INFO_fields_and_subfields['SVTYPE'] != 'BND'):
					# if POS and END are the same or differ by 1 for any SVTYPE except BND, 
					# and if REF has nucleotides, then use them to determine the true END
					if (int(END) <= (int(POS) + 1)): 
						if ( (int(POS) + len(REF) - 1) > int(END) ):
							END = str(int(POS) + len(REF) - 1 + 5)
							outline_part1 = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL

			# The end field still may not have the true end, eg. in the case of INDELs.
			# If there is more information to derive the end field, then use it.

			end_1 = int(END)
			end_2 = 0
			end_3 = 0
			end_4 = 0
			#if ('SVLEN' in INFO_fields_and_subfields):
			#	this_SVLEN = abs(int(INFO_fields_and_subfields['SVLEN'])) # SVLEN can be +ve or -ve
			#	end_2 = int(POS) + this_SVLEN - 1
			if ('SVEND' in INFO_fields_and_subfields):
				this_SVEND = abs(int(INFO_fields_and_subfields['SVEND']))
				end_3 = int(this_SVLEN)
			if ('END' in INFO_fields_and_subfields):
				this_END = abs(int(INFO_fields_and_subfields['END']))
				end_4 = int(this_END)
			max_end = max(end_1, end_2, end_3, end_4)
			if (max_end > end_1):
				END = str(max_end)
				outline_part1 = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL

			# If this is a mobile element insertion, then it may have had END filled in above using code meant for structural variants.
			# Reset END for mobile element insertions here, and set it to either POS or to the length of the TSD (tandem site duplication).

			if (len(all_ALTs) > 8):
				if (all_ALTs[0:8] == '<INS:ME:'):
					END = str(int(POS) + 1)
					outline_part1 = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL
					if ('TSD' in INFO_fields_and_subfields):
						if (INFO_fields_and_subfields['TSD'] != 'null'): # MELT produces null
							tsd_length = len(INFO_fields_and_subfields['TSD'])
							END = str(int(POS) + tsd_length)
							outline_part1 = CHROM + "\t" + POS + "\t" + END + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL

			# INFO_field_count_of_subfields[INFO_header]

			for i in range( 0, len(INFO_field_order) ): # For every info field in header, write out this variant's value for it
				this_info_key = INFO_field_order[i]
				this_info_value = ''
				if this_info_key in INFO_fields_and_subfields:
					this_info_value = str(INFO_fields_and_subfields[this_info_key])
				outline_part2 = outline_part2 + "\t" + this_info_value

			if (len(infields) > 9) : # For every sample, write out its first field, which is probably the GT field
				for i in range( 9, len(infields) ):
					this_sample = str(infields[i])
					this_sample_fields = this_sample.split(":")
					this_sample_GT = str(this_sample_fields[0])
					outline_part2 = outline_part2 + "\t" + this_sample_GT

			this_VAF_in_INFO = '.' # write out the INFO VAF field if there is one
			if ('VAF' in INFO_fields_and_subfields):
				this_VAF_in_INFO = str(INFO_fields_and_subfields['VAF'])
			outline_part2 = outline_part2 + "\t" + this_VAF_in_INFO

			if (len(infields) > 9) : # For every sample, write out its sample_format VAF field if there is one
				# If this is manta, then there is no VAF but we can derive it from PR and SR in the FORMAT fields:
				# ##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
				# ##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999">
				# GT:FT:GQ:PL:PR:SR	0/1:PASS:387:750,0,384:21,13:17,14
				# GT:FT:GQ:PL:PR:SR	1/1:MinGQ:4:139,5,0:0,0:0,6
				# GT:FT:GQ:PL:PR:SR       0/0:HomRef:143:0,93,602:21,0:34,0       0/0:HomRef:155:0,105,674:26,0:35,0      0/1:PASS:225:275,0,405:16,6:26,8
				# GT:FT:GQ:PL:PR:SR       0/1:PASS:70:474,0,67:0,0:8,15   0/1:PASS:113:647,0,110:0,0:11,17        0/1:PASS:77:772,0,74:0,0:11,23
				# GT:FT:GQ:PL:PR:SR       1/1:PASS:88:999,91,0:0,0:0,40   1/1:PASS:108:999,111,0:0,0:0,51 1/1:PASS:106:999,109,0:0,0:0,44
				# GT:FT:GQ:PL:PR:SR       1/1:PASS:33:662,36,0:0,0:0,32   0/1:PASS:75:158,0,72:2,0:6,10   0/1:PASS:73:147,0,70:1,0:6,10
				# However, this could be a non-manta input and it could have PR and SR that mean something else, so verify before using PR and SR.
				vaf_idx = -1
				PR_idx = -1
				SR_idx = -1
				this_format = str(infields[8])
				this_format_fields = this_format.split(":")
				for i in range( 0, len(this_format_fields) ):
					this_format = this_format_fields[i]
					if (this_format == "VAF"):
						vaf_idx = i
					if (this_format == "PR"):
						PR_idx = i
					if (this_format == "SR"):
						SR_idx = i
				for i in range( 9, len(infields) ):
					this_sample_VAF = '.'
					this_sample = str(infields[i])
					this_sample_fields = this_sample.split(":")
					if (vaf_idx > -1):
						this_sample_VAF = str(this_sample_fields[vaf_idx])
					else:
						if ((PR_idx > -1) and (SR_idx > -1)):
							this_PR = str(this_sample_fields[PR_idx])
							this_SR = str(this_sample_fields[SR_idx])
							this_PR_bits = this_PR.split(",")
							this_SR_bits = this_SR.split(",")
							if ( (len(this_PR_bits) == 2) and (len(this_SR_bits) == 2) ):
								if ( (is_integer(this_PR_bits[0])) and (is_integer(this_PR_bits[1])) and (is_integer(this_SR_bits[0])) and (is_integer(this_SR_bits[1])) ):
									vaf_denominator = int(this_PR_bits[0]) + int(this_PR_bits[1]) + int(this_SR_bits[0]) + int(this_SR_bits[1])
									if (vaf_denominator > 1):
										vaf_enumerator = int(this_PR_bits[1]) + int(this_SR_bits[1])
										this_sample_VAF = float(vaf_enumerator) / float(vaf_denominator)
										this_sample_VAF = "%0.8f" % this_sample_VAF
										this_sample_VAF = str(this_sample_VAF)
					outline_part2 = outline_part2 + "\t" + this_sample_VAF

			outline = outline_part1 + outline_part2
			out_TSV.write(outline + "\n")

	out_TSV.close()

if __name__=='__main__':
    main()

