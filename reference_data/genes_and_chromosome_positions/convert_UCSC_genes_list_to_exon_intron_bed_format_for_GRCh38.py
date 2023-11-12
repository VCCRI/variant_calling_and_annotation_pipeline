#!/usr/bin/python
# python3 split_UCSC_refseq_genes_to_genes_exons_introns.py -i input_genes_file -out_genes output_genes_file -out_exons output_exons_file -out_introns out_introns_file -out_cds output_CDS_file -out_cdsexons output_CDSexons_file -out_utr5 output_utr5_file -out_utr3 output_utr3_file

# This program reads in UCSC genes, and outputs separate files containing genes, exons, and introns.
# Input and output files are in 1-based bed format.

# head UCSC_GRCh38_GenesAndGenePredictions_transcriptExonsOnOneLine_RefSeq_20200219.txt
# #chrom  strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        name2
# chr1    -       67092164        67134970        67093579        67127240        9       67092164,67096251,67103237,67111576,67115351,67125751,67127165,67131141,67134929,    67093604,67096321,67103382,67111644,67115464,67125909,67127257,67131227,67134970,    C1orf141
# chr1    -       67092164        67134970        67093004        67127240        8       67092164,67095234,67096251,67115351,67125751,67127165,67131141,67134929,     67093604,67095421,67096321,67115464,67125909,67127257,67131227,67134970,        C1orf141

# python3 ../programs/convert_UCSC_genes_list_to_exon_intron_bed_format_for_GRCh38.py \
#   -i UCSC_GRCh38_GenesAndGenePredictions_transcriptExonsOnOneLine_RefSeq_20200219.txt \
#   -out_genes UCSC_GRCh38_GenesAndGenePredictions_genes_RefSeq_20200219.txt \
#   -out_exons UCSC_GRCh38_GenesAndGenePredictions_exons_RefSeq_20200219.txt \
#   -out_introns UCSC_GRCh38_GenesAndGenePredictions_introns_RefSeq_20200219.txt \
#   -out_cds UCSC_GRCh38_GenesAndGenePredictions_CDS_RefSeq_20200219.txt \
#   -out_cdsexons UCSC_GRCh38_GenesAndGenePredictions_CDSexons_RefSeq_20200219.txt \
#   -out_utr5 UCSC_GRCh38_GenesAndGenePredictions_UTR5_RefSeq_20200219.txt \
#   -out_utr3 UCSC_GRCh38_GenesAndGenePredictions_UTR3_RefSeq_20200219.txt

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2020, Victor Chang Cardiac Research Institite'
# This program is a modification of a program by the same author at Garvan Institute of Medical Research and Kinghorn Cancer Centre

import sys
import os
import subprocess
import argparse

######################################################
def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

######################################################
def main():

	parser = argparse.ArgumentParser(description='Read in UCSC table text file of genes, output bed files for genes, exons and introns.')
	parser.add_argument('-i', action="store", dest="in_UCSC", required=True, help='Input UCSC table text file list of genes')
	parser.add_argument('-out_genes', action="store", dest="out_genes", required=True, help='Output bed format file of genes')
	parser.add_argument('-out_exons', action="store", dest="out_exons", required=True, help='Output bed format file of exons')
	parser.add_argument('-out_introns', action="store", dest="out_introns", required=True, help='Output bed format file of introns')
	parser.add_argument('-out_cds', action="store", dest="out_cds", required=True, help='Output bed format file of CDS')
	parser.add_argument('-out_cdsexons', action="store", dest="out_cdsexons", required=True, help='Output bed format file of CDS exons')
	parser.add_argument('-out_utr5', action="store", dest="out_utr5", required=True, help='Output bed format file of UTR-5')
	parser.add_argument('-out_utr3', action="store", dest="out_utr3", required=True, help='Output bed format file of UTR-3')
	args = parser.parse_args()

	in_UCSC = open(args.in_UCSC, 'r')
	inlines = in_UCSC.readlines()
	in_UCSC.close()

	out_genes = open(args.out_genes, 'w')
	out_exons = open(args.out_exons, 'w')
	out_introns = open(args.out_introns, 'w')
	out_cds = open(args.out_cds, 'w')
	out_cdsexons = open(args.out_cdsexons, 'w')
	out_utr5 = open(args.out_utr5, 'w')
	out_utr3 = open(args.out_utr3, 'w')

	for i in range( 1, len(inlines) ): # skip the first line, it is the header line
		inline = inlines[i]
		stripped_inline = inline.strip()
		if (stripped_inline != ''):
			infields = inline.split("\t")
			infields[ len(infields)-1 ] = infields[ len(infields)-1 ].strip() # remove the carriage return on last field
			in_chrom = str(infields[0])
			#if (len(input_field_chrom) >= 3):
			#	if (input_field_chrom[0:3] == 'chr'):
			#		input_field_chrom = input_field_chrom[3:]
			in_strand = str(infields[1])
			in_txStart = str(infields[2])
			in_txEnd = str(infields[3])
			in_cdsStart = str(infields[4])
			in_cdsEnd = str(infields[5])
			in_exonCount = str(infields[6])
			in_exonStarts = str(infields[7])
			in_exonEnds = str(infields[8])
			in_geneSymbol = str(infields[9])

			outline_gene = in_chrom + "\t" + in_txStart + "\t" + in_txEnd + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
			out_genes.write( outline_gene )
			outline_cds = in_chrom + "\t" + in_cdsStart + "\t" + in_cdsEnd + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
			out_cds.write( outline_cds )
			if (int(in_txStart) < int(in_cdsStart)):
				outline_utr5 = in_chrom + "\t" + in_txStart + "\t" + in_cdsStart + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
				out_utr5.write( outline_utr5 )
			if (int(in_txEnd) > int(in_cdsEnd)):
				outline_utr3 = in_chrom + "\t" + in_cdsEnd + "\t" + in_txEnd + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
				out_utr3.write( outline_utr3 )

			in_exon_starts = in_exonStarts.split(",")
			in_exon_ends = in_exonEnds.split(",")
			in_cdsStart = int(in_cdsStart)
			in_cdsEnd = int(in_cdsEnd)
			for j in range( 0, len(in_exon_starts) ):
				this_exon_start = in_exon_starts[j]
				this_exon_end = in_exon_ends[j]
				if (this_exon_start != ''):
					this_exon_start = int(this_exon_start)
					this_exon_end = int(this_exon_end)
					outline_exons = in_chrom + "\t" + str(this_exon_start) + "\t" + str(this_exon_end) + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
					out_exons.write( outline_exons )
					this_exon_falls_inside_cds = False
					this_cdsexon_start = this_exon_start
					this_cdsexon_end = this_exon_end
					if ((in_cdsStart <= this_exon_start) and (in_cdsEnd >= this_exon_end)):
						this_exon_falls_inside_cds = True
					else:
						if ((this_exon_start <= in_cdsStart) and (in_cdsStart <= this_exon_end)):
							this_exons_falls_inside_cds = True
							this_cdsexon_start = in_cdsStart
						if ((this_exon_end >= in_cdsEnd) and (in_cdsEnd >= this_exon_start)):
							this_exon_falls_inside_cds = True
							this_cdsexon_end = in_cdsEnd
					if (this_exon_falls_inside_cds):
						outline_cdsexons = in_chrom + "\t" + str(this_cdsexon_start) + "\t" + str(this_cdsexon_end) + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
						out_cdsexons.write( outline_cdsexons )

			first_exon_starts = str(in_exon_starts[0])
			last_exon_ends = str(in_exon_ends[ len(in_exon_ends)-1 ])
			if (last_exon_ends == ''):
				last_exon_ends = str(in_exon_ends[ len(in_exon_ends)-2 ])
			if ((in_txStart != first_exon_starts) or (in_txEnd != last_exon_ends)):
				print( 'Unexpected mismatch of transcript and exon starts and ends for', in_txStart, first_exon_starts, in_txEnd, last_exon_ends, 'inline', inline )

			for j in range( 1, len(in_exon_starts) ):
				if (in_exon_starts[j] != ''):
					this_intron_start = str( int(in_exon_ends[j-1]) + 1 )
					this_intron_end = str( int(in_exon_starts[j]) - 1 )
					if (int(this_intron_start) < int(this_intron_end)):
						outline_introns = in_chrom + "\t" + this_intron_start + "\t" + this_intron_end + "\t" + in_geneSymbol + "\t" + in_strand + "\n"
						out_introns.write( outline_introns )

	out_genes.close()
	out_exons.close()
	out_introns.close()
	out_cds.close()
	out_cdsexons.close()
	out_utr5.close()
	out_utr3.close()

if __name__=='__main__':
    main()


