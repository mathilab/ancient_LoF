#Created 10 June 2021 by Nina Tansey
#Description: List LoF variants in Neanderthal exomes
#Parameters: argv[1] = gnomad.v2.1.1.all_lofs.txt.bgz 
#  argv[2:] = chromosome VCF nonref files 
#    e.g. /content/Archaics_chr22_mq25_mapab100_nonref_only.vcf.gz

import re
import sys
import gzip

def LoF_dictionary(LoF_path): 
  header = 0
  chrom_index = -1 
  pos_index = -1 
  ref_index = -1
  alt_index = -1
  gene_ids_index = -1
  gene_symbols_index = -1
  LoF_dict = {}

  #build the LoF dataframe from the gnomad dataset
  with gzip.open(LoF_path, 'rt') as f:
    for line in f:
      values = line.split()
      if header == 0: #header line, set indexes 
        for i in range(0, len(values)):
          if values[i] == 'chrom': chrom_index = i
          if values[i] == 'pos': pos_index = i
          if values[i] == 'ref': ref_index = i
          if values[i] == 'alt': alt_index = i
          if values[i] == 'gene_ids': gene_ids_index = i
          if values[i] == 'gene_symbols': gene_symbols_index = i
        #check indexes set
        if -1 in [chrom_index, pos_index, ref_index, alt_index, gene_ids_index, 
                  gene_symbols_index]: 
          print('error: no header line to locate gene_symbols')
        header = 1
      #add variant to dict: (chrom, pos) = (ref, alt, gene_ids, gene_symbols)
      LoF_dict[(str(values[chrom_index]), values[pos_index])] = (
               values[ref_index], values[alt_index], values[gene_ids_index],
                values[gene_symbols_index])
  return LoF_dict

def print_variants(chrom_path, LoF_dict):
  chrom_index = 0 #standard VCF format
  pos_index = 1 #standard VCF format
  ref_index = 3 #standard VCF format
  alt_index = 4 #standard VCF format
  
  with gzip.open(chrom_path,'rt') as f:
      for line in f:
        #check for comment line (##)
        if line[:2] == '##': continue 
        #header
        if line[:1] == '#': continue 
        values = line.split()
        #variant present in dictionary
        key = (values[chrom_index], values[pos_index])
        if key in LoF_dict.keys():
          #check the reference and alternate alleles match
          if (values[ref_index] == LoF_dict[key][0]) and (values[alt_index] == LoF_dict[key][1]):            
            for i in range(9, len(values)):
              #comment this if block to print only heterozygous LoF
              if values[i] == '1/1':
                #uncomment this to print transitions
                # if ((values[ref_index] in ['A', 'G'] and values[alt_index] in ['A', 'G']) or (
                #   values[ref_index] in ['C', 'T'] and values[alt_index] in ['C', 'T'])):
                #uncomment this to print transversions
                # if not ((values[ref_index] in ['A', 'G'] and values[alt_index] in ['A', 'G']) or (
                #   values[ref_index] in ['C', 'T'] and values[alt_index] in ['C', 'T'])):
                print(values[pos_index], '\t', values[chrom_index], '\t', 
                  values[ref_index], '\t', values[alt_index], '\t', 
                  LoF_dict[key][2], '\t', LoF_dict[key][3]) 
                break
              #comment this if block to print only homozygous LoF
              if values[i] in ['0/1', '1/0']:
                #uncomment this to print transitions
                # if ((values[ref_index] in ['A', 'G'] and values[alt_index] in ['A', 'G']) or (
                #   values[ref_index] in ['C', 'T'] and values[alt_index] in ['C', 'T'])):
                #uncomment this to print transversions
                # if not ((values[ref_index] in ['A', 'G'] and values[alt_index] in ['A', 'G']) or (
                #   values[ref_index] in ['C', 'T'] and values[alt_index] in ['C', 'T'])):
                print(values[pos_index], '\t', values[chrom_index], '\t', 
                 values[ref_index], '\t', values[alt_index], '\t', 
                 LoF_dict[key][2], '\t', LoF_dict[key][3])
                break
  return 

def main( ):
  #paths for VCF files
  args = sys.argv[2:]

  LoF_path = sys.argv[1]
  LoF_dict = LoF_dictionary(LoF_path)

  #print header line for txt file
  print('pos\tchrom\tref\talt\tgene_ids\tgene_symbols')
  
  for path in args:
    print_variants(path, LoF_dict)

if __name__ == "__main__":
    main()