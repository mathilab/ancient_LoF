#Created 03 June 2021 by Nina Tansey
#Description: Count LOF homozygous/heterozygous variants in Neanderthal xsomes
#Parameters: argv[1] = gnomad.v2.1.1.all_lofs.txt.bgz 
#  argv[2:] = chromosome VCF nonref files 
#    e.g. /content/Archaics_chr22_mq25_mapab100_nonref_only.vcf.gz

import pandas as pd
import re
import sys

def LoF_dictionary(LoF_path): 
  #build the LoF dataframe from the gnomad dataset
  LoF_df = pd.read_csv(LoF_path, compression='gzip', header=0, sep='\t', 
                      quotechar='"')
  #build the LoFs variant dictionary from the LoF dataframe
  LoFs = {}
  for i in range (0, LoF_df.shape[0]): 
    #variant = {(chrom,pos): (ref,alt,gene_ids,gene_symbols)}
    LoFs[(str(LoF_df['chrom'][i]), LoF_df['pos'][i])] = (LoF_df['ref'][i], 
      LoF_df['alt'][i], LoF_df['gene_ids'][i], LoF_df['gene_symbols'][i], 
      LoF_df['chrom'][i])
  return LoFs

def count_variants(chrom_path, LoFs, count):
  #make dataframe for Neanderthal data (this is just chromosome 22)
  #make header dynamic instead of hard coded TODO, there's a comment tag I can use
  #can also read file line by line
  chrom_df = pd.read_csv(chrom_path, compression='gzip', header=12, sep='\t', 
    quotechar='"')

  #initialize variant counts
  Altai = count[0]
  Vindija = count[1]
  Denisova = count[2]

  #find the chromosome number in the filepath
  temp = re.findall(r'\d+', chrom_path)
  chrom_num = str(list(map(int, temp))[0])

  #print(chrom_num)

  # look at the variants present
  for i in range(0, chrom_df.shape[0]): #while not EOF -- can use read_line() : something like that
      if (chrom_num,chrom_df['POS'][i]) in LoFs.keys():
        #print('POS= ', chrom_df['POS'][i])
        if (chrom_df['REF'][i] is LoFs[(chrom_num,chrom_df['POS'][i])][0][0]) and (
            chrom_df['ALT'][i] is LoFs[(chrom_num,chrom_df['POS'][i])][1][0]):
          #print('variant at REF=' + chrom_df['REF'][i] + ', ALT=' + chrom_df['ALT'][i])
          print(LoFs[(chrom_num,chrom_df['POS'][i])])
          #print('POS= ', chrom_df['POS'][i])
          #print(chrom_df['AltaiNeandertal'][i], ', ', chrom_df['Vindija33.19'][i], ', ', chrom_df['Denisova'][i])
          #check where the variants are
          #print('same REF[0] and ALT[0]') # don't do this, want lengths of each to be 1
          #print(' ')
          #count Altai
          if (chrom_df['AltaiNeandertal'][i] in ['0/1', '1/0']):
            Altai['heterozygous'] += 1
          if (chrom_df['AltaiNeandertal'][i] == '1/1'):
            Altai['homozygous'] += 1
          #count Vindija
          if (chrom_df['Vindija33.19'][i] in ['0/1', '1/0']):
            Vindija['heterozygous'] += 1
          if (chrom_df['Vindija33.19'][i] == '1/1'):
            Vindija['homozygous'] += 1
          #count Denisova
          if (chrom_df['Denisova'][i] in ['0/1', '1/0']):
            Denisova['heterozygous'] += 1
          if (chrom_df['Denisova'][i] == '1/1'):
            Denisova['homozygous'] += 1
          #print('')
  #Return the variant counts per Neanderthal genome
  #print('Altai: ', Altai)
  #print('Vindija: ', Vindija)
  #print('Denisova: ', Denisova)
  return [Altai, Vindija, Denisova]

def main( ):
  #paths for chromosome files
  #'/content/Archaics_chr22_mq25_mapab100_nonref_only.vcf.gz'
  args = sys.argv[2:]

  #print('count')
  LoF_path = sys.argv[1]
  LoFs = LoF_dictionary(LoF_path)

  #initialize count
  Altai = {'homozygous': 0, 'heterozygous': 0}
  Vindija = {'homozygous': 0, 'heterozygous': 0}
  Denisova = {'homozygous': 0, 'heterozygous': 0}
  count = [Altai, Vindija, Denisova]
  
  for path in args:
    count = count_variants(path, LoFs, count)
  print('Altai: ', count[0])
  print('Vindija: ', count[1])
  print('Denisova: ', count[2])

if __name__ == "__main__":
    main()