#!/usr/bin/env python3
mtgrasp_version = 'v1.1.7'

import argparse
import sys
import shlex
import subprocess
import os
import pandas as pd
from Bio import SeqIO
import numpy as np
import re

parser = argparse.ArgumentParser(description='Usage of mtgrasp_standardize.py')
parser.add_argument('-i', '--input', help='Input fasta file containing your mitochondrial sequence', required=True)
parser.add_argument('-o', '--output', help='Output directory', required=True)
parser.add_argument('-p', '--prefix', help='Prefix for the output file', default='mtgrasp_standardized')
parser.add_argument('-c', '--gencode', help='Mitochondrial genetic code used for gene annotation', required=True)
parser.add_argument('-a', '--annotate', help='Run gene annotation on the final assembly output [False]', action='store_true')
parser.add_argument('-mp', '--mitos_path', help='Complete path to runmitos.py', default = None)

args = parser.parse_args()
file = args.input
mito_gencode = args.gencode
output_dir = args.output
sample = args.prefix
annotate = args.annotate
mitos_path = args.mitos_path


# get the directory of the script
string = subprocess.check_output(['which', 'mtgrasp_standardize.py'])
string = string.decode('utf-8')
# remove new line character
string = string.strip()
# split string by '/'
script_dir = '/'.join(string.split('/')[0:-1])


# This function checks if the tRNA-Phe gene is on the forward or reverse strand
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

# This functions checks the strand orientation of the assembly sequence when the tRNA-Phe gene is not present
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

    
    if match_list.count(True) >= match_list.count(False):
       return 'forward'
    else:
       return 'reverse'

   

# This function checks if the tRNA-Phe gene is present in the assembly
def check_if_trnF_gaa_in_fasta(fas_file):
    id_list = []
    for record in SeqIO.parse(fas_file, "fasta"):
        id_list.append(str(record.description))
    
    if '(gaa)' in ''.join(id_list):
        return True
    else:
        return False

# This function finds the start and end position of the tRNA-Phe gene
def find_trnF_gaa_pos(fas_file):
    id_list = []
    for record in SeqIO.parse(fas_file, "fasta"):
        id_list.append(str(record.description))
    header = ''.join([v for v in id_list if '(gaa)'in v])
    start_pos = int(''.join(header.split(';')[1].split(' ')[1].split('-')[0]))
    end_pos = int(''.join(header.split(';')[1].split(' ')[1].split('-')[1]))
    return start_pos, end_pos


# This function standardizes the start site of the assembly sequence if it is on the forward strand
def forward_strand_start_site_standardization(seq, start_pos):
    return seq[start_pos-1:] + seq[0:start_pos-1]

# This function standardizes the start site of the assembly sequence if it is on the reverse strand
def reverse_strand_start_site_standardization(seq, end_pos):
    return seq[end_pos:] + seq[0:end_pos]

# This function reverse complements the assembly sequence
def reverse_complement(seq):
    return seq.translate(str.maketrans('ATCGatcg', 'TAGCtagc'))[::-1]

# This function removes the duplicates in a list, 
# this is to remove duplicated sequences after start-site standardization,
#  when both scenario 1 and scenario 2 sequences consist of non-fragmented tRNA-Phe gene
def remove_duplicates_in_a_list(list):
    # convert the list to an array
    list = np.array(list)
    # drop the duplicates
    list = np.unique(list)
    # convert the array to a list
    return list.tolist()

# This function writes the standardized sequence to a fasta file
def write_fasta_file(file_name, standardized_seq, standardization_status, sample):
  
    if len(standardized_seq) > 1:
      for i, seq_item in enumerate(standardized_seq): 
        with open(file_name, 'a') as f:
          f.write('>%s_%s_seq_%s\n' % (standardization_status, sample, i+1) + seq_item + '\n')
    else:
      with open(file_name, 'a') as f:
        f.write('>%s_%s_seq\n' % (standardization_status, sample) + standardized_seq[0] + '\n')


# This function finds the path to the conda environment
def find_conda_env(env_name):
    # Find the environment path using conda info
    output = subprocess.run(["conda", "info", "--envs"], stdout=subprocess.PIPE, text=True).stdout
    pattern = fr"\s+(\S+envs/{env_name})$"
    match = re.search(pattern, output, re.MULTILINE)
    if match:
        print(f"Conda environment {env_name} found.")
        return match.group(1)
    else:
        # print the error message and exit
        print(f"Conda environment {env_name} not found.")
        sys.exit(1)


        

# This function runs MitoS annotation
def run_mitos(env_name, file, code, dir, script_dir, mitos_path):
    # Run mitos via conda
    if mitos_path == 'None':
       path_to_env = find_conda_env(env_name)
       cmd = f"conda run -p {path_to_env} runmitos.py -i {file} --noplots  -c {code} -o {dir} --linear --refdir {script_dir}/data/refseqs_mitos -r refseq81m"
       process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
       output = process.communicate()[0].decode("utf-8").strip()
       return output
    # Run mitos independently without conda
    else:
      if os.path.exists(f"{mitos_path}/runmitos.py"):

        # Add mitos_path to PATH
        new_path = mitos_path
        current_path = os.environ.get('PATH', '')
        new_path_value = f'{new_path}:{current_path}'
        os.environ['PATH'] = new_path_value

        cmd = f"{mitos_path}/runmitos.py -i {file} --noplots  -c {code} -o {dir} --linear --refdir {script_dir}/data/refseqs_mitos -r refseq81m"
        process = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
        output = process.communicate()[0].decode("utf-8").strip()
        return output
      else:
        return "The provided mitos path is not correct, please double check"
    



print('Start standardization for %s' % (file))
    
# make a new directory for the annotation of each assembly
anno_dir = f'{output_dir}/annotations'
cmd1 = f'mkdir -p {anno_dir}' 
cmd_shlex1 = shlex.split(cmd1)
subprocess.call(cmd_shlex1)

#check if file empty
num_of_assemblies = len([record for record in SeqIO.parse(file, "fasta")])
if 'Scenario' in open(file).read():
    scenario_in_file = True
else:
    scenario_in_file = False
# if the file is empty, skip the start-site standardization
if os.stat(file).st_size == 0:
    print('File is empty, no start-site standardization needed')
    # write an empty fasta file
    with open('%s/%s.final-mtgrasp_%s-assembly.fa'%(output_dir,sample, mtgrasp_version), 'w') as f:
         f.write('')


# elif the file has multiple sequences and 'Scenario' not in file, 
# skip the annotation, because the start-site standardization is not needed in this case
# start-site standardization is only needed when there is only one sequence in the fasta file
elif num_of_assemblies > 1 and scenario_in_file == False:
    print('Multiple contigs in the fasta file, standardization cannot be performed')
    # write the original fasta file to the output directory
    with open('%s/%s.final-mtgrasp_%s-assembly.fa'%(output_dir,sample, mtgrasp_version), 'w') as f:
            for record in SeqIO.parse(file, "fasta"):
                f.write('>Non-Standardized_%s_seq\n' % (sample) + str(record.seq) + '\n')
# if the file has only one sequence or has 2 sequences (but 'Scenario' is found in file)
# start the start-site standardization
elif num_of_assemblies == 2 or num_of_assemblies == 1: 
      print('Start standardizing the start site of the mitochondrial genome')
      code = mito_gencode
      print('The mitochondrial genetic code is: %s' % (code))
        
      
      seq_list = []
      for record in SeqIO.parse(file, "fasta"):
        seq_list.append(str(record.seq)) 
      num_seq = len(seq_list)

      # annotate the file
      try:
          output = run_mitos("mitos", file, code, anno_dir, script_dir, mitos_path)
          print(f"Output: {output}")
      except Exception as e:
          print(f"Error: {e}")
      
      standardized_seq = []
      if num_seq > 1 and 'Scenario' not in open(file).read():
         print("'Scenario' not found in file when expected, please check the input fasta file")  # This is because 'Scenario' string should be found in the one-piece assembly file that has been gap-filled
        
      elif num_seq > 1 and 'Scenario' in open(file).read():
            trnF_seq_length = 0
            for i in range(num_seq):
                fas_file = '%s/%s/result.fas'%(anno_dir, i)
         
                if check_if_trnF_gaa_in_fasta(fas_file):
                   print('tRNA-Phe Gene Found!')
                   start_pos, end_pos = find_trnF_gaa_pos(fas_file)
                   if end_pos - start_pos > trnF_seq_length:
                      trnF_seq_length = end_pos - start_pos
                      seq_index = i
                      start_trna = start_pos
                      end_trna = end_pos
                else:
                    print('tRNA-Phe Gene Not Found!')
                    if non_trnF_strand_check(fas_file) == 'forward':
                       standardized_seq.append(seq_list[0])
                    else:
                       standardized_seq.append(reverse_complement(seq_list[0]))
            if trnF_seq_length != 0:
               if trnF_strand_check(fas_file) == 'forward':
                   standardized_seq.append(forward_strand_start_site_standardization(seq_list[seq_index], start_trna))
               else:
                   standardized_seq.append(reverse_complement(reverse_strand_start_site_standardization(seq_list[seq_index], end_trna)))

      
      else:
          fas_file = '%s/result.fas'%(anno_dir)
          if check_if_trnF_gaa_in_fasta(fas_file) == True:
              print('tRNA-Phe Gene Found!')
              start_pos, end_pos = find_trnF_gaa_pos(fas_file)
              if trnF_strand_check(fas_file) == 'forward':
                  standardized_seq.append(forward_strand_start_site_standardization(seq_list[0], start_pos))
              else:
                  standardized_seq.append(reverse_complement(reverse_strand_start_site_standardization(seq_list[0], end_pos)))
              
              
          else:
              print('tRNA-Phe Gene Not Found!')
              if non_trnF_strand_check(fas_file) == True:
                   standardized_seq.append(seq_list[0])
              else:
                   standardized_seq.append(reverse_complement(seq_list[0]))

              
      
      standardized_seq = remove_duplicates_in_a_list(standardized_seq)
      file_name = '%s/%s.final-mtgrasp_%s-assembly.fa'%(output_dir, sample, mtgrasp_version)

      
      # Store the results of the function and file reading in variables
      trnF_gaa_present = check_if_trnF_gaa_in_fasta(fas_file)
      file_contents = open(file).read()
      if 'Scenario' in file_contents:
        scenario = True
      else:
        scenario = False

      if trnF_gaa_present and scenario:
           write_fasta_file(file_name, standardized_seq, 'StartSite_Strand_Standardized_Circular', sample)
      elif trnF_gaa_present and not scenario:
           write_fasta_file(file_name, standardized_seq, 'StartSite_Strand_Standardized_Linear', sample)
      elif not trnF_gaa_present and scenario:
           write_fasta_file(file_name, standardized_seq, 'Strand_Standardized_Circular', sample)
      else:
           write_fasta_file(file_name, standardized_seq, 'Strand_Standardized_Linear', sample)
else:
        print('Error: Something wrong with the fasta file')
        sys.exit(1)
    
# Remove MitoS annotation intermediate files
cmd1 = f'rm -r {output_dir}/annotations'
cmd_shlex1 = shlex.split(cmd1)
subprocess.call(cmd_shlex1)

# Run MitoS annotation on the final assembly fasta file if '-a' argument is provided
if annotate:
      print('Start Annotating %s/%s.final-mtgrasp_%s-assembly.fa'%(output_dir, sample, mtgrasp_version))
      if not os.path.exists(f'{output_dir}/annotation_output'):
         cmd1 = f'mkdir -p {output_dir}/annotation_output' 
         cmd_shlex1 = shlex.split(cmd1)
         subprocess.call(cmd_shlex1)
      
      
      
      try:
          output = run_mitos("mitos", '%s/%s.final-mtgrasp_%s-assembly.fa'%(output_dir, sample, mtgrasp_version), mito_gencode,f'{output_dir}/annotation_output', script_dir, mitos_path)
          print(f"Output: {output}")
      except Exception as e:
          print(f"Error: {e}")
