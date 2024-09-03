#!/usr/bin/env python3

'''
This script can be used to summarize mtGrasp assembly outputs by providing a text file containing the relative or complete path(s) to assembly output folder(s). 
'''
mtgrasp_version = 'v1.1.5'

import argparse
import os
import os.path
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help='Input text file containing the path(s) to assembly output folder(s)',required=True)
parser.add_argument('-p','--prefix',help='Prefix for output files',required=True)

args = parser.parse_args()

input = args.input
prefix = args.prefix

tsv_filename = f'{prefix}_mtgrasp_{mtgrasp_version}_assembly_summary.tsv'
path_txt_filename = f'{prefix}_mtgrasp_{mtgrasp_version}_path_to_output.txt'


# Count the number of contigs,  Ns per 1000bp, and number of bases, length of longest contig, standardization status, etc. in an assembly output fasta
def get_assembly_metrics(fasta):
  seq_count = 0
  num_bases = 0
  total_n_count = 0 #Ns per 1000bp
  with open(fasta, 'r') as f:
      seq_len_list = []
      for line in f:
         line = line.strip()
         if line.startswith('>'):
              fasta_header = line
         else:
              seq = line.upper()
              seq_len = len(seq)
              seq_len_list.append(seq_len)
              num_bases += seq_len
              total_n_count += seq.count("N") # number of Ns per fasta sequence
              seq_count += 1
         
      if num_bases != 0:
           n_count_per_1000 = total_n_count /num_bases/ 1000
      else:
           n_count_per_1000 = 0
      max_seq_len = max(seq_len_list)
      if seq_count == 0:
          standardization_status = 'N/A'
          circle_check = 'N/A'
      elif seq_count == 1:
          if "Circular" in fasta_header:
              circle_check = 'Circular'
          else:
              circle_check = 'Linear'

          if "StartSite" in fasta_header and "Strand" in fasta_header:
              standardization_status = "StartSite_Strand_Standardized"
          elif "StartSite" in fasta_header and "Strand" not in fasta_header:
            standardization_status = "StartSite_Standardized"
          elif "StartSite" not in fasta_header and "Strand" in fasta_header:
              standardization_status = "Strand_Standardized"
          else:
              standardization_status = "Non-Standardized"
      else:
          # Circularization and standardization are only conducted for one-piece assemblies to avoid potentially introducing misassemblies 
          circle_check = 'Linear'
          standardization_status = "Non-Standardized"


      return n_count_per_1000, seq_count, num_bases, max_seq_len, circle_check, standardization_status 
      
assembly_list, n_count_list, seq_count_list, num_bases_list, max_seq_list, circle_check_list, standardization_status_list, file_list = [], [], [], [], [], [], [], []

with open(input, 'r') as dirs:
    for outdir in dirs:
        outdir = outdir.strip()
        if not os.path.exists(outdir):
            print(f"Error: mtGrasp output directory {outdir} does not exist")
        directory = "%s/final_output"%(outdir)
        dir_exists=os.path.exists(directory)

        if dir_exists:
            for dir in os.listdir(directory):

                fasta = os.path.abspath("%s/%s/%s.final-mtgrasp_%s-assembly.fa"%(directory,dir,dir, mtgrasp_version)) # Complete path to output fasta file
                # check if fasta output file exists 
                if os.path.exists(fasta) and any(line.startswith('>') for line in open(fasta)):
                    n_count_list.append(get_assembly_metrics(fasta)[0])
                    seq_count_list.append(get_assembly_metrics(fasta)[1])
                    num_bases_list.append(get_assembly_metrics(fasta)[2])
                    max_seq_list.append(get_assembly_metrics(fasta)[3])
                    circle_check_list.append(get_assembly_metrics(fasta)[4])
                    standardization_status_list.append(get_assembly_metrics(fasta)[5])
                else:
                    n_count_list.append(0)
                    seq_count_list.append(0)
                    num_bases_list.append(0)
                    max_seq_list.append(0)
                    circle_check_list.append('N/A')
                    standardization_status_list.append('N/A')
                
                
                assembly_list.append(dir)
                file_list.append(fasta)


df = pd.DataFrame({"Assembly": assembly_list,
                   "Ns per 1000bp": n_count_list,
                   "Number of Contigs": seq_count_list,
                   "Total Number of Base Pairs Per Assembly": num_bases_list,
                   "Length of the Longest Contig (bp)": max_seq_list,
                   "Circular or Linear": circle_check_list,
                   "Standardization Status": standardization_status_list
                   })

df.to_csv(tsv_filename, sep="\t", index=False)

with open(path_txt_filename, 'w') as f:
    for path in file_list:
        f.write(f"{path}\n")


             
