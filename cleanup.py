#!/usr/bin/env python3

# This python script can be used to clean up unneeded intermediate files to free up your disk space


import sys
import shlex
import subprocess


dir = sys.argv[1] # the output directory user specified for mtGasp (-o or --out_dir)


cmd = f'cd {dir} && rm -r -v !(final_*|*.gz)' 
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)