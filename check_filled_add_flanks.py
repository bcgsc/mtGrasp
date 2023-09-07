#!/usr/bin/env python


'''
This script checks if the fake gap is filled by sealer. If it is, the script adds the gap-filled construct to the original assembly.

'''


import sys
import shlex
import subprocess
import os
from Bio import SeqIO


sealer_log = sys.argv[1]
fake_gap_filled = sys.argv[2]
assembly = sys.argv[3]
output = sys.argv[4]
pre_sealer_fake_gap = sys.argv[5]


with open(sealer_log) as file_in:
            # save the lines in the sealer log file to a list
            ls = []
            for l in file_in:
                ls.append(l)
# check if the gap was filled
gap_remaining = int(float(ls[1].split(' ')[0]) - float(ls[-3].split(' ')[3])) # gap remaining after sealer
pre_sealer_gap = int(float(ls[1].split(' ')[0])) # gap before sealer
if pre_sealer_gap - gap_remaining == 1:
# check the length of the fasta sequence
    for record in SeqIO.parse(pre_sealer_fake_gap, "fasta"):
        # calculate the flank size
        flank_size = int((len(str(record.seq))-10)/2)

        cmd = f"add_flanks.py {fake_gap_filled} {assembly} {flank_size} {output}"
        args = shlex.split(cmd)
        subprocess.call(args)
# if the gap was not filled, copy the original assembly to the output
else:
    cmd = f"cp {assembly} {output}"
    args = shlex.split(cmd)
    subprocess.call(args)
