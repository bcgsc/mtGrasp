#!/usr/bin/env python3

# Split a fasta file into multiple files, one per sequence and create a reference config file for ntjoin
# Usage: split_fasta.py <fasta file> <output directory> <config file>

import sys

ntjoin_reference_weight = 2 #Refer to ntJoin README for more info: https://github.com/bcgsc/ntJoin
input = sys.argv[1]
output = sys.argv[2]
config = sys.argv[3]
ref_list = []
with open(input, 'r') as f:
    for line in f:
        if line[i].startswith('>'):
            header = lines[i].strip()
            filename = header.split(' ')[0].strip('>') + '.fasta'
            sequence = lines[i+1].strip()
            with open(output + '/' + filename, 'w') as g:
                g.write(header + '\n')
                g.write(sequence + '\n')
            ref_list.append(f'{output}/{filename},{ntjoin_reference_weight}')

with open(config, 'w') as f:
    f.write('\n'.join(ref_list))
                

