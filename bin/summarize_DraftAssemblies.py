#!/usr/bin/env python3

# This script generates a summary of the assembly outputs (including the number of one-pieces, the maximum scaffold length in each assembly, and the assembly name) 
# and writes the list of path to assembly output files to a text file.



import argparse
import os
import os.path
from Bio import SeqIO
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',help='Input file that contains the path to each library folder',required=True)
parser.add_argument('-c','--output',help='Output csv file',required=True)
parser.add_argument('-l','--pathlist',help='Output a list of paths to the output files',required=True)

args = parser.parse_args()

input = args.input
output = args.output
pathlist = args.pathlist



def summarize_outputs(input, output, pathlist):
    assembly = []
    number_of_onepieces = []
    scaffold_max_lengths = []
    file_list =[]
    gap_remaining = []
    pre_sealer_gap = []
    sealer_prefix = []

    file = open(input, 'r')
    lines = file.read().splitlines()
    

    for line in lines:
        #summarize assembly outputs
        directory = "%s/assembly_output"%(line)
        dir_exists=os.path.isdir(directory)
        if dir_exists == True:
            for file in os.listdir(directory):
                file_list.append(os.path.abspath("%s/%s"%(directory,file)))
                basename = os.path.basename(file)
                assembly.append(basename.split('.')[0])
                number_of_onepieces.append(len([1 for line in open("%s/%s"%(directory,file)) if line.startswith(">")]))
                max_len = 0

                for record in SeqIO.parse("%s/%s"%(directory,file), "fasta"):
                    if len(record) > max_len:
                        max_len = len(record)
                scaffold_max_lengths.append(max_len)
    
    
        #summarize sealer log files
        dir_sealer = "%s/sealer"%(line)
        dir_sealer_exists=os.path.isdir(dir_sealer)
        if dir_sealer_exists == True:
            for f in os.listdir(dir_sealer):
                f_basename = os.path.basename(f)
                if f_basename.endswith("_log.txt") and os.stat("%s/%s"%(dir_sealer,f)).st_size != 0:
                    sealer_prefix.append("%s_%s"%(os.path.basename(line),f_basename.split('.')[0]))
                    with open("%s/%s"%(dir_sealer,f)) as file_in:
                       ls = []
                       for l in file_in:
                         ls.append(l)
                    gap_remaining.append(int(float(ls[1].split(' ')[0]) - float(ls[-3].split(' ')[3])))
                    pre_sealer_gap.append(int(float(ls[1].split(' ')[0])))

    df_assm = pd.DataFrame({'Assembly':assembly, 'Number of one-pieces':number_of_onepieces, 'Maximum scaffold length':scaffold_max_lengths})
    df_sealer = pd.DataFrame({'Assembly':sealer_prefix, 'Number of gaps (pre-GapFilling)':pre_sealer_gap, 'Number of gaps (post-GapFilling)':gap_remaining})
                        
    df = pd.merge(df_assm, df_sealer, on='Assembly')


    
    df.to_csv(output)
    with open(pathlist, 'w') as f:
            for item in file_list:
                f.write("%s\n" % item)
            
                


if __name__ == '__main__':
    print(summarize_outputs(input, output, pathlist))
