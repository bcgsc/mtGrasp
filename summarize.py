#!/usr/bin/env python3

# This script generates a summary of the assembly outputs (including the number of one-pieces, the maximum scaffold length in each assembly, and the assembly name) 
# and writes the list of path to assembly output files to a text file.



import argparse
import os
import os.path
from Bio import SeqIO
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help='Input text file containing the path(s) to assembly output folder(s)',required=True)
parser.add_argument('-t','--output',help='Output tsv filename',required=True)
parser.add_argument('-l','--pathlist',help='Output file consisting of the full path(s) to the assembly output file(s)',required=True)

args = parser.parse_args()

input = args.input
output = args.output
pathlist = args.pathlist



def summarize_outputs(input, output, pathlist):
    # report the sequence lengths of each scaffold in a fasta file
    def count_seq_lengths(file):
       lengths = []
       with open(file, "r") as f:
         for line in f:
            if line.startswith(">"):
                continue
            else:
                lengths.append(len(line.strip()))
       return ",".join([str(x) for x in lengths])

# report the standardization status of each assembly
    def extract_headers(file):
      headers = []
      with open(file, "r") as f:
        for line in f:
            if line.startswith(">"):
                line = line.strip()[1:] #remove >
                line = line.split("_")[:-4] #exclude values after 4th last underscore
                headers.append("_".join(line))
      return ",".join(headers)


    assembly = []
    number_of_sequences = []
    scaffold_lengths = []
    file_list =[]
    gap_remaining = []
    pre_sealer_gap = []
    sealer_prefix = []
    status = []

    file = open(input, 'r')
    lines = file.read().splitlines()
    

    for line in lines:
        #summarize sealer log files
        dir_sealer = "%s/prepolishing"%(line)
        dir_sealer_exists=os.path.isdir(dir_sealer)
        if dir_sealer_exists == True:
            for f in os.listdir(dir_sealer):
                f_basename = os.path.basename(f)
                if f_basename.endswith("_log.txt") and os.stat("%s/%s"%(dir_sealer,f)).st_size != 0:
                    sealer_prefix.append("%s_%s"%(os.path.basename(line),'_'.join(f_basename.split('.')[0].split('_')[:-1])))
                    with open("%s/%s"%(dir_sealer,f)) as file_in:
                       ls = []
                       for l in file_in:
                         ls.append(l)
                    gap_remaining.append(int(float(ls[1].split(' ')[0]) - float(ls[-3].split(' ')[3])))
                    pre_sealer_gap.append(int(float(ls[1].split(' ')[0])))
        else:
            print("Error: %s does not exist"%(dir_sealer))
            continue

        
        #summarize assembly outputs
        directory = "%s/standardized_output"%(line)
        dir_exists=os.path.exists(directory)
        if dir_exists == True:
            for dir in os.listdir(directory):
                fasta = os.path.abspath("%s/%s/post_standardization.fasta"%(directory,dir))
                file_list.append(fasta)
                
                assembly.append(dir)
                n = len([1 for l in open(fasta) if l.startswith(">")])
                number_of_sequences.append(n)
                scaffold_lengths.append(count_seq_lengths(fasta))
                status.append(extract_headers(fasta))
        else:
            print("Error: %s does not exist"%(directory))
            continue
    
    
        
    df_assm = pd.DataFrame({'Assembly':assembly, 'Number of sequences':number_of_sequences, 'Scaffold lengths':scaffold_lengths, 'Standardization status':status})
    df_sealer = pd.DataFrame({'Assembly':sealer_prefix, 'Number of gaps (pre-GapFilling)':pre_sealer_gap, 'Number of gaps (post-GapFilling)':gap_remaining})
                        
    df = pd.merge(df_assm, df_sealer, on='Assembly')


    
    df.to_csv(output,sep='\t',index=False)
    with open(pathlist, 'w') as f:
            for item in file_list:
                f.write("%s\n" % item)
            
                


if __name__ == '__main__':
    print(summarize_outputs(input, output, pathlist))
