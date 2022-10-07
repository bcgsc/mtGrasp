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
            ls = []
            for l in file_in:
                ls.append(l)
gap_remaining = int(float(ls[1].split(' ')[0]) - float(ls[-3].split(' ')[3]))
pre_sealer_gap = int(float(ls[1].split(' ')[0]))
if pre_sealer_gap - gap_remaining == 1:
# check the length of the fasta sequence
    for record in SeqIO.parse(pre_sealer_fake_gap, "fasta"):
        flank_size = int((len(str(record.seq))-10)/2)
        cmd = f"python3 /projects/transabyss/btl/itrack_dna/assemblies/bin/add_flanks.py {fake_gap_filled} {assembly} {flank_size} {output}"
        args = shlex.split(cmd)
        subprocess.call(args)
else:
    cmd = f"cp {assembly} {output}"
    args = shlex.split(cmd)
    subprocess.call(args)
