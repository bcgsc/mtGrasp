#!/usr/bin/env python3

# Split a fasta file into multiple files, one per sequence and create a reference config file for ntJoin
# Usage: split_fasta.py <fasta file> <output directory> <config file>

import sys
from itertools import islice

ntjoin_reference_weight = 2 #Refer to ntJoin README for more info: https://github.com/bcgsc/ntJoin
input = sys.argv[1]
output = sys.argv[2]
config = sys.argv[3]
ref_list = []
with open(input, 'r') as f:
    for i, line in enumerate(f):
        if line.startswith('>'):
            header = line.strip()
            filename = header.split(' ')[0].strip('>') + '.fasta'
            sequence = next(islice(f, 1, None), None).strip()
            with open(output + '/' + filename, 'w') as g:
                g.write(header + '\n')
                g.write(sequence + '\n')
            ref_list.append(f'{output}/{filename},{ntjoin_reference_weight}')

with open(config, 'w') as f:
    f.write('\n'.join(ref_list))
                

