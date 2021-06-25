#Created: 25 June 2021 by Nina Tansey
#Description: Build a dataset of synonymous and LoF variants and their 
# frequenciesfrom LoF list, where each syn variant is in a LoF gene
#Parameters: ARGV[1] = lof txt file, ARGV[2] = gnomad vcf dataset
#   e.g. gnomad.exomes.r2.1.1.sites.vcf.bgz

import gzip 
import argparse

def make_lofs_dict(lofs_path):
  LoF_dict = {} #LoF_dict[gene_id] = [(chrom,pos,ref,alt)]
  genes_dict = {} #genes_dict[gene_id] = chrom
  with open(lofs_path,'r') as f:
      for line in f:
        values = line[:-1].split()
        #account for header
        if values[0] == "CHROM" or values[0] == "chrom" or "#" in values[0]: continue
        if values[0] == "POS" or values[0] == "pos" or "#" in values[0]: continue
        #add LoF entry to dict
        if values[4] in LoF_dict.keys():
          LoF_dict[values[4]].append((values[1], values[0], values[2], values[3]))
        else:
          LoF_dict[values[4]] = [(values[1], values[0], values[2], values[3])]
        #add gene_id entry to dict
        genes_dict[values[4]] = values[1]
  return (LoF_dict, genes_dict)

def print_freqs(data_path, lofs_dict, genes_dict):
  #really should be dynamically retrieving indexes :(
  #syns_dict[gene_id] = ([logAlleleFrequency,chrom,pos,ref,alt,dif])
  with gzip.open(data_path, 'rt') as f:
    for line in f: 
      values = line[:-1].split()
      #account for header
      if "#" in values[0]: continue
      for gene_id in genes_dict.keys():
        #look at genes only in relevant chromosome
        if genes_dict[gene_id] == values[0] and gene_id in values[7]: 
          freq = values[7].split(";")[2]
          if "synonymous_variant" in values[7]: 
            print("%s\t%s\t%s\t%s\t%s\t%s\tsynonymous"%(values[0], values[1], 
                values[3], values[4], gene_id, freq))
          for i in range(0, len(lofs_dict[gene_id])):
            pos = lofs_dict[gene_id][i][1]
            ref = lofs_dict[gene_id][i][2]
            alt = lofs_dict[gene_id][i][3]
            if pos == values[1] and ref == values[3] and alt == values[4]:
              print("%s\t%s\t%s\t%s\t%s\t%s\tLoF"%(values[0], pos, ref, alt, 
                gene_id, freq))
  return 

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("lof_file", help="lof_filepath")
  parser.add_argument("data_file", help="data_filepath")                 
  args = parser.parse_args()

  dicts = make_lofs_dict(args.lof_file)

  print("#CHROM\tPOS\tREF\tALT\tGENE_ID\tALLELE_FREQ\tVARIANT")
  print_freqs(args.data_file, dicts[0], dicts[1])

if __name__ == "__main__":
    main()
