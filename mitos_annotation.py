#!/usr/bin/env python3

import sys
import shlex
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import numpy as np
import re


file = sys.argv[1] # input fasta file
mito_gencode = sys.argv[2] # mitochondrial genetic code used for gene annotation
output_dir = sys.argv[3] # output directory for mitos annotation
sample = sys.argv[4] # unique identifier for the assembly sample (e.g., accession number, species name or library ID number etc.)


# get the directory of the script
string = subprocess.check_output(['which', 'mitos_annotation.py'])
string = string.decode('utf-8')
# remove new line character
string = string.strip()
# split string by '/'
script_dir = '/'.join(string.split('/')[0:-1])
   
def trnF_strand_check(fas_file):
  id_list = []
  for record in SeqIO.parse(fas_file, "fasta"):
          id_list.append(str(record.description))
  header = ''.join([v for v in id_list if '(gaa)'in v])
  if '; +;' in header:
    return 'forward'
  else:
    return 'reverse'

def check_if_match(x, y):
    if x == y:
        return True
    else:
        return False

def non_trnF_strand_check(fas_file):
    

    mito_genes_df = pd.read_csv(f'{script_dir}/data/mito_gene_orientation.csv')
    id_list = []
    for record in SeqIO.parse(fas_file, "fasta"):
        id_list.append(str(record.description))
    orientation = []
    gene_name = []

    for id in id_list:
       gene_name.append(''.join(id.split(';')[-1].split(' ')[1]))
       orientation.append(''.join(id.split(';')[-2].split(' ')[1]))

    gene_orientation = pd.DataFrame({'gene': gene_name, 'orientation': orientation})

    # merge the two dataframes
    mito_genes_df = pd.merge(mito_genes_df, gene_orientation, on='gene', how='left')

    # drop rows with nan values
    mito_genes_df = mito_genes_df.dropna().reset_index(drop=True)
    mito_genes_df['orientation'] = mito_genes_df['orientation'].apply(lambda x: 1 if x == '+' else -1)
    mito_genes_df['match'] = mito_genes_df.apply(lambda x: check_if_match(x['orientation'], x['strand']), axis=1)
    match_list = mito_genes_df['match'].tolist()

    

    if match_list.count(True) > match_list.count(False):
       return 'forward'
    elif match_list.count(True) == match_list.count(False):
       return 'forward'
    else:
       return 'reverse'

   


def check_if_trnF_gaa_in_fasta(fas_file):
    id_list = []
    for record in SeqIO.parse(fas_file, "fasta"):
        id_list.append(str(record.description))
    
    if '(gaa)' in ''.join(id_list):
        return True
    else:
        return False

def find_trnF_gaa_pos(fas_file):
    id_list = []
    for record in SeqIO.parse(fas_file, "fasta"):
        id_list.append(str(record.description))
    header = ''.join([v for v in id_list if '(gaa)'in v])
    start_pos = int(''.join(header.split(';')[1].split(' ')[1].split('-')[0]))
    end_pos = int(''.join(header.split(';')[1].split(' ')[1].split('-')[1]))
    return start_pos, end_pos



def forward_strand_start_site_standardization(seq, start_pos):
    return seq[start_pos-1:] + seq[0:start_pos-1]
def reverse_strand_start_site_standardization(seq, end_pos):
    return seq[end_pos:] + seq[0:end_pos]
def reverse_complement(seq):
    return seq.translate(str.maketrans('ATCGatcg', 'TAGCtagc'))[::-1]

def remove_duplicates_in_a_list(list):
    # convert the list to an array
    list = np.array(list)
    # drop the duplicates
    list = np.unique(list)
    # convert the array to a list
    return list.tolist()

def write_fasta_file(file_name, standardized_seq, standardization_status, sample):
  
    if len(standardized_seq) > 1:
      for i in range(len(standardized_seq)):
        with open(file_name, 'a') as f:
          f.write('>%s_%s_seq_%s\n' % (standardization_status, sample, i+1) + standardized_seq[i] + '\n')
    else:
      with open(file_name, 'a') as f:
        f.write('>%s_%s_seq\n' % (standardization_status, sample) + standardized_seq[0] + '\n')



def find_conda_env(env_name):
    # Find the environment path using conda info
    output = subprocess.run(["conda", "info", "--envs"], stdout=subprocess.PIPE, text=True).stdout
    pattern = fr"\s+(\S+envs/{env_name})$"
    match = re.search(pattern, output, re.MULTILINE)
    if match:
        print(f"Conda environment {env_name} found.")
        return match.group(1)
    else:
        print(f"Conda environment {env_name} not found. Installing Mitos to a conda environment called mitos")
        subprocess.run(["conda", "install", "mitos", "-c", "bioconda", "-m", "-n", "mitos"])
        find_conda_env(env_name)



def run_mitos(env_name, file, code, anno_dir, script_dir):
    path_to_env = find_conda_env(env_name)

    cmd = f"conda run -p {path_to_env} runmitos.py -i {file} -c {code} -o {anno_dir} --linear --refdir {script_dir}/data/refseqs_mitos -r refseq81m"
    process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    output = process.communicate()[0].decode("utf-8").strip()

    return output



print('Annotating file: %s' % (file))
    
# make a new directory for the annotation of each assembly
anno_dir = output_dir
cmd1 = f'mkdir -p {anno_dir}' 
cmd_shlex1 = shlex.split(cmd1)
subprocess.call(cmd_shlex1)

#check if file empty
# if the file is empty, skip the annotation
if os.stat(file).st_size == 0:
    print('File is empty, no annotation')
    # write an empty fasta file
    with open('%s/post_standardization.fasta'%(output_dir), 'w') as f:
         f.write('')


    
# if file is not empty, annotate
else: 
      
      code = mito_gencode
      print('The mitochondrial genetic code is: %s' % (code))
        
      # get the number of sequences in the fasta file
      num_seq = len(list(SeqIO.parse(file, "fasta"))) 
      seq_list = []
      for record in SeqIO.parse(file, "fasta"):
        seq_list.append(str(record.seq))     
      # annotate the file
      
      shell_name = os.environ['SHELL'].split('/')[-1]
      
      try:
          output = run_mitos("mitos", file, code, anno_dir, script_dir)
          print(f"Output: {output}")
      except Exception as e:
          print(f"Error: {e}")
      
      standardized_seq = []
      if num_seq > 1:
         for i in range(num_seq):
            fas_file = '%s/%s/result.fas'%(anno_dir, i)

            if check_if_trnF_gaa_in_fasta(fas_file) == True:
              start_pos, end_pos = find_trnF_gaa_pos(fas_file)
              if trnF_strand_check(fas_file) == 'forward':
                  standardized_seq.append(forward_strand_start_site_standardization(seq_list[i], start_pos))
              else:
                  standardized_seq.append(reverse_complement(reverse_strand_start_site_standardization(seq_list[i], end_pos)))
        
            else:
              if non_trnF_strand_check(fas_file) == 'forward':
                   standardized_seq.append(seq_list[i])
              else:
                   standardized_seq.append(reverse_complement(seq_list[i]))
      
      else:
          fas_file = '%s/result.fas'%(anno_dir)
          if check_if_trnF_gaa_in_fasta(fas_file) == True:
              start_pos, end_pos = find_trnF_gaa_pos(fas_file)
              if trnF_strand_check(fas_file) == 'forward':
                  standardized_seq.append(forward_strand_start_site_standardization(seq_list[0], start_pos))
              else:
                  standardized_seq.append(reverse_complement(reverse_strand_start_site_standardization(seq_list[0], end_pos)))
              
              
          else:
              if non_trnF_strand_check(fas_file) == True:
                   standardized_seq.append(seq_list[0])
              else:
                   standardized_seq.append(reverse_complement(seq_list[0]))
              
      
      standardized_seq = remove_duplicates_in_a_list(standardized_seq)
      file_name = '%s/post_standardization.fasta'%(output_dir)

      
      if check_if_trnF_gaa_in_fasta(fas_file) == True and 'Scenario' in open(file).read():
           write_fasta_file(file_name, standardized_seq, 'StartSite_Strand_Standardized_EndRecovered', sample)
      elif check_if_trnF_gaa_in_fasta(fas_file) == True and 'Scenario' not in open(file).read():
           write_fasta_file(file_name, standardized_seq, 'StartSite_Strand_Standardized_NotEndRecovered', sample)
      elif check_if_trnF_gaa_in_fasta(fas_file) == False and 'Scenario' in open(file).read():
           write_fasta_file(file_name, standardized_seq, 'Strand_Standardized_EndRecovered', sample)
      else:
           write_fasta_file(file_name, standardized_seq, 'Strand_Standardized_NotEndRecovered', sample)
