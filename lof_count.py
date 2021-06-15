#Created 03 June 2021 by Nina Tansey
#Description: Count LOF homozygous/heterozygous variants in Neanderthal exomes
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

def count_variants(chrom_path, LoF_dict, samples):
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
        if key in LoF_dict.keys():
          #check the reference and alternate alleles match
          if (values[ref_index] == LoF_dict[key][0]) and (values[alt_index] == LoF_dict[key][1]):
            for i in range(9, len(values)):
              if values[i] == '1/1':
                samples[i-9]['homozygous'] += 1
              if values[i] in ['0/1', '1/0']:
                samples[i-9]['heterozygous'] += 1
  return samples

def main( ):
  #paths for VCF files
  args = sys.argv[2:]

  LoF_path = sys.argv[1]
  LoF_dict = LoF_dictionary(LoF_path)

  #initialize counts
  counts = []
  
  for path in args:
    counts = count_variants(path, LoF_dict, counts)
  for i in range(0, len(counts)):
    print(counts[i])

if __name__ == "__main__":
    main()