#Created 03 June 2021 by Nina Tansey
#Description: Count syn homozygous/heterozygous variants in Neanderthal exomes
#Parameters: argv[1] = gnomad.v2.1.1.all_lofs.txt.bgz 
#  argv[2:] = chromosome VCF nonref files 
#    e.g. /content/Archaics_chr22_mq25_mapab100_nonref_only.vcf.gz

import re
import sys
import gzip
import argparse

def syn_dictionary(syn_path): 
  header = 0
  chrom_index = -1 
  pos_index = -1 
  ref_index = -1
  alt_index = -1
  gene_ids_index = -1
  syn_dict = {}

  #build the syn dataframe from the gnomad dataset
  with open(syn_path, 'rt') as f:
    for line in f:
      values = line.split()
      if header == 0: #header line, set indexes 
        #might have issues setting chrom index :(
        for i in range(0, len(values)):
          if values[i] == 'chrom' or values[i] == '#CHROM': chrom_index = i
          if values[i] == 'pos' or values[i] == 'POS': pos_index = i
          if values[i] == 'ref' or values[i] == 'REF': ref_index = i
          if values[i] == 'alt' or values[i] == 'ALT': alt_index = i
          if values[i] == 'gene_ids' or values[i] == "GENE_ID": gene_ids_index = i
        #check indexes set
        if -1 in [chrom_index, pos_index, ref_index, alt_index, gene_ids_index]: 
          print('error: no header line to locate gene_symbols')
        header = 1
      #add variant to dict: (chrom, pos) = (ref, alt, gene_ids)
      syn_dict[(str(values[chrom_index]), values[pos_index])] = (
               values[ref_index], values[alt_index], values[gene_ids_index])
  return syn_dict

def count_variants(chrom_path, syn_dict, samples):
  chrom_index = 0 #standard VCF format
  pos_index = 1 #standard VCF format
  ref_index = 3 #standard VCF format
  alt_index = 4 #standard VCF format
  
  with gzip.open(chrom_path,'rt') as f:
      for line in f:
        #check for comment line (##)
        if line[:2] == '##': continue 
        #header
        if line[:1] == '#': 
          #samples = list of dicts {(index): (sample, homozygous, heterozygous)}
          if not samples:
            values = line[:-1].split()
            for i in range(9, len(values)):
              samples.append({'name': values[i], 'homozygous': 0, 'heterozygous': 0})
          continue
        values = line.split()
        #variant present in dictionary
        key = (values[chrom_index], values[pos_index])
        if key in syn_dict.keys():
          #check the reference and alternate alleles match
          if (values[ref_index] == syn_dict[key][0]) and (values[alt_index] == syn_dict[key][1]):
            for i in range(9, len(values)):
              if values[i] == '1/1':
                samples[i-9]['homozygous'] += 1
              if values[i] in ['0/1', '1/0']:
                samples[i-9]['heterozygous'] += 1
  return samples

def main( ):
  args = sys.argv[2:]

  syn_path = sys.argv[1]
  syn_dict = syn_dictionary(syn_path)

  #initialize counts
  counts = []
  
  for path in args:
    counts = count_variants(path, syn_dict, counts)
  for i in range(0, len(counts)):
    print(counts[i])

if __name__ == "__main__":
    main()