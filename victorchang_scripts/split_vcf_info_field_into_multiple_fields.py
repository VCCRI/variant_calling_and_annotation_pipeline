#!/usr/bin/python
# python3 split_vcf_info_field_into_multiple_fields.py -i input_VCF_file -o input_VCF_file -f field_name -d delimiter

# python3 split_vcf_info_field_into_multiple_fields.py -i test_19W001217_multinucleotide.spliceai_tensorflow.vcf -o test_19W001217_multinucleotide.spliceai_tensorflow.spliceai_info.vcf -f SpliceAI -d '|'

# ##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant annotation. These include delta scores (DS) and delta positions (DP) for acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">

# VARIANT_TYPE=intronic_._;ALLELE_END;SpliceAI=T|NEXN|0.00|0.00|0.00|0.00|-39|3|-6|-39	GT:AD:DP:GQ:PL	0/1:15,15:30:99:586,0,2087


__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2020, Victor Chang Medical Research Institute'

import sys
import os
import argparse
import re

######################################################
def is_integer(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

######################################################
def split_info_data( old_info, field_name, delimiter, info_subfield_names ):

	old_info = str(old_info)
	new_info = old_info
	new_subfields = ''
	if (old_info != '.'):
		info_fields = old_info.split(';')
		for this_info_field in info_fields:
			bits = this_info_field.split('=')
			if (len(bits) >= 2):
				this_key = str(bits[0])
				this_value = str(bits[1])
				if (this_key == field_name):
					info_subfield_values = this_value.split(delimiter)
					for i in range( 0, len(info_subfield_names) ):
						this_subfield_name = str(info_subfield_names[i])
						this_subfield_value = str(info_subfield_values[i])
						new_subfield_pair = this_subfield_name + '=' + this_subfield_value
						if (new_subfields == ''):
							new_subfields = new_subfield_pair
						else:
							new_subfields = new_subfields + ';' + new_subfield_pair
		new_info = new_info + ';' + new_subfields

	return new_info

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in an uncompressed VCF file having INFO fields annotations, and output this record with INFO fields as tab-separated columns in annovar data format, with first five columns: Chr Start End Ref Alt.')
	parser.add_argument('-i', action="store", dest="in_VCF", required=True, help='Input VCF file')
	parser.add_argument('-o', action="store", dest="out_VCF", required=True, help='Output VCF file')
	parser.add_argument('-f', action="store", dest="field_name", required=True, help='This is the info field that needs to be split')
	parser.add_argument('-d', action="store", dest="delimiter", required=False, help='This is the delimiter that the field should be split on')
	args = parser.parse_args()

	# read each record of input VCF file and process it

	in_VCF = open(args.in_VCF, 'r')
	out_VCF = open(args.out_VCF, 'w')
	field_name = 'SpliceAI'
	if (args.field_name is not None):
		field_name = str(args.field_name)
	delimiter = '|'
	if (args.delimiter is not None):
		delimiter = str(args.delimiter)

	info_subfield_names = []

	in_header = True

	for inline in in_VCF:

		inline = inline.strip()

		if (in_header == True):

			# process each header line read in

			if (inline[0:1] == '#'):
				in_header = True
				if (len(inline) > 11):

					# Process the INFO header lines.

					if (inline[0:11] == '##INFO=<ID='):

						info_regex = re.compile('##INFO=<ID=(.*)>')
						this_info_definition = info_regex.search(inline).group(1)
						bits_by_comma = this_info_definition.split(',')
						this_info_name = str(bits_by_comma[0])

						if (this_info_name == field_name):

							info_regex2 = re.compile('Format: (.*)">')
							this_info_definition2 = info_regex2.search(inline).group(1)
							this_info_sub_names = this_info_definition2

							info_subfields = this_info_sub_names.split(delimiter)
							for this_subfield in info_subfields:
								new_subfield_name = this_info_name + '_' + str(this_subfield)
								info_subfield_names.append( new_subfield_name )

					elif (inline[0:6] == '#CHROM'):

						for this_subfield in info_subfield_names:
							new_line = '##INFO=<ID=' + this_subfield + ',Number=.,Type=String,Description="' + field_name + ' subfield">'
							out_VCF.write( new_line + "\n" )

				out_VCF.write( inline + "\n" )

			else:
				in_header = False

		if (in_header == False):

			infields = inline.split("\t")
			CHROM = str(infields[0])
			POS = str(infields[1])
			SNV_ID = str(infields[2])
			REF = str(infields[3])
			all_ALTs = str(infields[4])
			QUAL = str(infields[5])
			FILTER = str(infields[6])
			INFO = str(infields[7])

			new_info = split_info_data( INFO, field_name, delimiter, info_subfield_names )

			outline = CHROM + "\t" + POS + "\t" + SNV_ID + "\t" + REF + "\t" + all_ALTs + "\t" + QUAL + "\t" + FILTER + "\t" + new_info

			if (len(infields) >= 9):
				for i in range( 8, len(infields) ):
					this_field = str(infields[i])
					outline = outline + "\t" + this_field

			out_VCF.write( outline + "\n" )

	out_VCF.close()


if __name__=='__main__':
    main()

