#Created: 16 June 2021 by Nina Tansey
#Description: Outputs LoF variants oe scores in a csv
#Parameters: argv[1] = .gz file of oe scores. e.g. gnomad.v2.1.1.lof_metrics.by_gene.txt.gz
    # argv[2] = Tabulated .txt file of LoF variants

import sys
import gzip

def make_LoFs_dict(LoF_path):
  count = 0
  index = -1
  LoFs_dict = {}
  with open(LoF_path,'r') as f:
      for line in f:
        #values stores the tabulated values from the line in a list
        values = line[:-1].split()
        #account for header
        if count == 0: 
          for i in range(0, len(values)):
            #print('i',i)
            #print(values[i])
            if values[i] == 'gene_symbols': index = i
          #print('index', index)
          if index == -1: print('error: no header line to locate gene_symbols')
          count = count + 1
          #print('count', count)
          continue
        #if count is 10: break
        #dynamically from first line get index for "gene_symbols" location
        LoFs_dict[values[index]] = ""
        #print('got line', line[:-1])
        
        #print(values[5])
        count = count + 1
  #print(LoFs_dict)
  return LoFs_dict

def update_oe(allLoFs_path, LoFs_dict):
  count = 0
  gene_index = -1
  oe_lof_index = -1
  all_LoFs = {}
  with gzip.open(allLoFs_path,'rt') as f:
      for line in f:
        #values stores the tabulated values from the line in a list
        values = line[:-1].split()
        #account for header
        if count == 0: 
          #dynamically find indexes
          for i in range(0, len(values)):
            if values[i] == 'oe_lof': oe_lof_index = i
            if values[i] == 'gene': gene_index = i
          if -1 in [oe_lof_index, gene_index]: print('error: no header line locate gene/oe_lof index')
          count = count + 1
          continue
        if values[gene_index] in LoFs_dict.keys():
          LoFs_dict[values[gene_index]] = values[oe_lof_index]
        all_LoFs[values[gene_index]] = values[oe_lof_index]
        count = count + 1
  return LoFs_dict

def write_oe_csv(filename, LoFs_dict):
  import csv
  with open(filename, 'w') as f:
      f.write("gene_symbols,oe_lof\n")
      for key in LoFs_dict.keys():
        if LoFs_dict[key] == '' or LoFs_dict[key] == 'NA': continue
        num = float(LoFs_dict[key])
        f.write("%s,%f\n"%(key,num))
      f.close

def main():
  #paths for LoF files
  args = sys.argv[2:]

  #paths for metric data
  allLoFs_path = sys.argv[1]
  
  for path in args:
    LoFs_dict = make_LoFs_dict(path)
    #write filename dynamically -- change .txt file to .csv
    filename = path[:-4] + '.csv'
    write_oe_csv(filename, update_oe(allLoFs_path, LoFs_dict))

if __name__ == "__main__":
    main()