# python3 annotate_vcf_with_mitimpact_http_calls.py -i invcf.vcf -o outvcf.vcf
# python3 annotate_vcf_with_mitimpact_http_calls.py -i 19W001352.chromAll.annovar_clinvar20200629.chromM.vcf -o 19W001352.chromAll.annovar_clinvar20200629.chromM.mitimpact.vcf

# This program reads in a vcf of variants on the mitochondrial chromosome.
# For each variant it calls MitImpact, and annotates the variant with the MitImpact results for this variant.

# MitImpact documentation is available at:
# https://mitimpact.css-mendel.it/description#annotations

# An example MitImpact call and results:
# result = requests.get('https://mitimpact.css-mendel.it/api/v2.0/search_allele/6253/T>A')

# Annotation sources for MitImpact
# https://mitimpact.css-mendel.it/description
# 
# Conservation
#  * PhyloP 100V
#  * PhastCons 100V
# Pathogenicity meta-predictors
#  * CAROL
#  * Condel
#  * COVEC WMV
#  * Meta-SNP
#  * MtoolBox
#  * APOGEE
#  * DEOGEN2 score
# Cancer-specific predictors
#  * PolyPhen2
#  * SIFT
#  * MutationAssessor
#  * CHASM
# Residue interaction
#  * EVmutation
#  * Site A-B InterP
#  * Site A-B IntraP
#  * DDG intra
#  * DDG intra interface
#  * DDG inter
# Pathogenicicty predictors
#  * PolyPhen2
#  * SIFT
#  * FatHmm
#  * FatHmmW
#  * PROVEAN
#  * MutationAssessor
#  * EFIN SP
#  * EFIN HD
#  * VEST
#  * PANTHER
#  * PhD-SNP
#  * MutationTaster
#  * CADD
#  * SNAP
#  * Mitoclass1
#  * SNPDryad
# Databases of Allelic frequency & Phenotype
#  * ClinVar ClinSig
#  * ClinVar ClnDBN
#  * ClinVar ClnAllele id
#  * ClinVar ClnDSDB
#  * MITOMAP allele
#  * MITOMAP phenotype
#  * MITOMAP hooplasmy
#  * MITOMAP heteroplasmy
#  * MITOMAP status
#  * MITOMAP NRef
#  * COSMIC 90 id
#  * dbSNP 153 id
# Pathogenic variants, CPD
#  * Frequency
#  * AA ref
#  * AA alt
#  * Aln pos
#  * RefSeq protein id
#  * Species name
#  * Ncbi taxon id


import sys
import csv
import subprocess
import argparse
import requests
import json


def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Annotate a vcf of mitochondrial variants with MitImpact annotations from http JSON calls.')
    parser.add_argument('-i', dest='invcf', 
                        help='Input VCF of mitochondrial variants')
    parser.add_argument('-o', dest='outvcf', 
                        help='Output VCF file annotated with MitImpact information')
    return parser.parse_args()


def read_invcf(input_file_path):
  with open(input_file_path) as invcf_handle:
    inlines = invcf_handle.readlines()
  inlines_without_trailing_carriage_return = []
  for inline in inlines:
    new_inline = inline[ 0 : len(inline)-1 ]
    inlines_without_trailing_carriage_return.append(new_inline)
  return inlines_without_trailing_carriage_return


def convert_key(in_key):
  in_key = str(in_key)
  out_key = in_key.replace("-", "_")
  out_key = "MitImpact_" + out_key
  return out_key


def obtain_mitimpact_headers():

  list_of_mitimpact_hdrs = []
  result = requests.get('https://mitimpact.css-mendel.it/api/v2.0/search_allele/6253/T>A')
  # Or maybe the most recent query syntax is: http://mitimpact.css-mendel.it/search_allele/6253/T>A
  if (result.status_code == 200):
    data = result.json()
    for key1 in data['variant']:
      list_of_mitimpact_hdrs.append(str(key1))

  return list_of_mitimpact_hdrs


def obtain_and_write_out_headers(invcf, outvcf):

  list_of_mitimpact_hdrs = []

  for inline in invcf:
    if (len(inline) >=1):
      if (inline[0:1] == '#'):

        is_last_header = False
        if (len(inline) >= 6):
          if (inline[0:6] == '#CHROM'):
            is_last_header = True

        if (is_last_header):

          # Get the mitimpact headers
          list_of_mitimpact_hdrs = obtain_mitimpact_headers()
          for key1 in list_of_mitimpact_hdrs:
            key2 = convert_key(key1)
            new_hdr = '##INFO=<ID=' + str(key2) + ',Number=.,Type=String,Description="' + str(key2) + ' annotation provided by MitImpact">'
            outline = str(new_hdr) + "\n"
            outvcf.write(outline)

          outline = str(inline) + "\n"
          outvcf.write(outline)

        else:

          outline = str(inline) + "\n"
          outvcf.write(outline)

  return list_of_mitimpact_hdrs


def obtain_mitimpact_annotation(chrom, pos, ref, alt, list_of_mitimpact_hdrs):

  dict_of_mitimpact_results = {}

  url = 'https://mitimpact.css-mendel.it/api/v2.0/search_allele/' + str(pos) + '/' + ref + '>' + alt
  result = requests.get(url)
  if (result.status_code == 200):

    data = result.json()

    for key1 in data['variant']:
      key_result = data['variant'][str(key1)]
      #key2 = convert_key(key1)
      dict_of_mitimpact_results[key1] = str(key_result)

  return dict_of_mitimpact_results


def add_info(old_info, list_of_mitimpact_hdrs, dict_of_mitimpact_results):

  new_info = old_info
  for key1 in list_of_mitimpact_hdrs:

    if (key1 in dict_of_mitimpact_results):
      key_result = dict_of_mitimpact_results[key1]
      key2 = convert_key(key1)
      key_output = str(key2) + '=' + str(key_result)

      if ((new_info == '.') or (new_info == '')):
        new_info = key_output
      else:
        new_info = new_info + ';' + key_output

  return new_info

def annotate_vcf_records(invcf, outvcf, list_of_mitimpact_hdrs):

  for inline in invcf:
    if (len(inline) >=1):
      char1 = inline[0:1]
      if (char1 != '#'):

        # process this input vcf data record
        infields = inline.split("\t")
        chrom = str(infields[0])
        pos = str(infields[1])
        ref = str(infields[3])
        alt = str(infields[4])
        old_info = infields[7]

        # obtain mitimpact annotation to put onto end of INFO field
        dict_of_mitimpact_results = obtain_mitimpact_annotation(chrom, pos, ref, alt, list_of_mitimpact_hdrs)

        final_info = add_info(old_info, list_of_mitimpact_hdrs, dict_of_mitimpact_results)
        infields[7] = final_info

        # output the vcf record with correctly formatted info field
        outline = str(infields[0])
        for i in range(1, len(infields)):
          outline = outline + "\t" + str(infields[i])
        outline = str(outline) + "\n"
        outvcf.write(outline)
  return


def convert_vep_multiple_transcript_vcf_to_official_vcf_format(input_file_path, output_file_path):

  invcf = read_invcf(input_file_path)

  outvcf = open(output_file_path, 'w')

  list_of_mitimpact_hdrs = obtain_and_write_out_headers(invcf, outvcf)

  annotate_vcf_records(invcf, outvcf, list_of_mitimpact_hdrs)

  return


def main(args):

  args = parse_arguments(args)
  input_file_path = args.invcf
  output_file_path = args.outvcf

  convert_vep_multiple_transcript_vcf_to_official_vcf_format(input_file_path, output_file_path)

  #print(' ')
  #print('End of convert_vep_multiple_transcript_vcf_to_official_vcf_format')


if __name__ == "__main__":
    main(sys.argv)


