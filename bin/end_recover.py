# Check if the assembly sequence has overlapping regions 
# Create fake gap construct using the flanking regions of the original assembly (e.g., end + NNNNNNNNNN + start)
# If the overlap is found, the length of the flanking regions is the length of overlaping region
# If there is no ovelap, the flanking region length will be set to 200bp

import sys
import shlex
import subprocess
import os

assembly = sys.argv[1]
bf = sys.argv[2]
p = sys.argv[3]


folder = ".".join(assembly.split("/")[-1].split(".")[0:-3])
dir = "/".join(assembly.split("/")[0:-2])
r1 = "/".join(assembly.split("/")[0:-2]) + "/" + "r1.fq.gz"
r2 = "/".join(assembly.split("/")[0:-2]) + "/" + "r2.fq.gz"

out_dir = dir + "/end_recovery/" + folder
cmd = f'mkdir -p {out_dir}'
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)


cmd = f'python3 bin/overlap_check.py {assembly} 1 {out_dir}/fake_gap_unfilled.fa'
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)


print("Start sealer gap filling")
cmd = f'abyss-sealer -b{bf} -k 40 -k 60 -k 80 -k 100 -P {p} -o {out_dir}/fake_gap_filled -S {out_dir}/fake_gap_unfilled.fa {r1} {r2}'
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)


print("Start adding flanks back to the original assembly") 
cmd = f'python3 bin/check_filled_add_flanks.py {out_dir}/fake_gap_filled_log.txt {out_dir}/fake_gap_filled_scaffold.fa {assembly} {out_dir}/flank_added_assembly.fa {out_dir}/fake_gap_unfilled.fa' 
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)

