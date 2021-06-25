#Created: 23 June 2021 by Nina Tansey
#Description: Build list of syn variants closest in freq to LoF vars from file 
#   of frequncies for synonymous and LoF variants
#Parameters: Dataset of variants and frequencies
#   e.g. arch_syn_freqs.gz = synonymous frequencies
#        arch_lof_freqs.gz = LoF frequencies

import gzip
import math
import argparse

def make_lofs_dict(lofs_path):
  LoF_freqs = {} #LoF_freqs[gene_id] = (logAlleleFrequency)
  with gzip.open(lofs_path,'rt') as f:
      for line in f:
        values = line[:-1].split()
        #account for header
        if "#" in values[0]: continue
        if values[4] in LoF_freqs.keys():
          LoF_freqs[values[4]].append(math.log(float(values[5][3:]),10))
        else:
          LoF_freqs[values[4]] = [math.log(float(values[5][3:]),10)]
  return LoF_freqs

def make_syns_dict(syns_path, lofs_freqs):
  #syns_dict[gene_id] = ([logAlleleFrequency,chrom,pos,ref,alt,dif])
  syns_dict = {} 
  with gzip.open(syns_path, 'rt') as f:
    for line in f: 
      values = line[:-1].split()
      #account for header
      if "#" in values[0]: continue
      #check you're actually getting the allele frequency field
      if values[5][0:3] != "AF=": continue
      #skip variant if the frequency is 0 or is LoF
      if float(values[5][3:]) == 0: continue
      if values[6] == 'LoF': continue
      freq = math.log(float(values[5][3:]),10)
      gene_id = values[4]
      #syn variant from an irrelevant gene
      if gene_id not in lofs_freqs.keys(): continue
      for i in range(0,len(lofs_freqs[gene_id])):
        diff = abs(lofs_freqs[gene_id][i] - freq)
        #add syn variant entry to the list for the corresponding gene_id
        if gene_id in syns_dict.keys():
          syns_dict[gene_id].append(((i,diff),values[0],values[1],values[2],values[3]))
        elif gene_id not in syns_dict.keys():
          syns_dict[gene_id] = [((i,diff),values[0],values[1],values[2],values[3])]
  #sort on the difference between lof freq and syn freq
  for gene_id in syns_dict.keys():
    syns_dict[gene_id] = sorted(syns_dict[gene_id], key=lambda x: x[0])
  return syns_dict

def write_syns_file(filename, lofs_dict, syns_dict, num):
  with open(filename, 'w') as f:
    #file header
    f.write("#CHROM\tPOS\tREF\tALT\tGENE_ID\n")
    for gene_id in syns_dict.keys():
      #print $num syn variants per lof variant
      for j in range(0,len(lofs_dict[gene_id])):
        count = 0
        for i in range(0,len(syns_dict[gene_id])):
          if i >= len(syns_dict[gene_id]):
            print("missing # of syns: ", num-i)
            break
          if syns_dict[gene_id][i][0][0] != j: continue
          if count == num: break
          f.write("%s\t%s\t%s\t%s\t%s\n"%(syns_dict[gene_id][i][1], 
                  syns_dict[gene_id][i][2], syns_dict[gene_id][i][3], 
                  syns_dict[gene_id][i][4], gene_id))
          count += 1
  

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-n", "--num_syn", type=int, 
                        help="increase output number of synonymous variants")
  parser.add_argument("lof_file", help="lof_filepath")
  parser.add_argument("syn_file", help="syn_filepath") 
  parser.add_argument("--out", help="output file name")                     
  args = parser.parse_args()

  if not args.lof_file or not args.syn_file: 
      print("error: missing file arguments")
      exit()
  lofs_dict = make_lofs_dict(args.lof_file)
  syns_dict = make_syns_dict(args.syn_file, lofs_dict)

  if args.num_syn:
      num = args.num_syn
  else: num=1

  if args.out: write_syns_file(args.out, lofs_dict, syns_dict, num)
  else: write_syns_file("output.txt", lofs_dict, syns_dict, num)

if __name__ == "__main__":
    main()
