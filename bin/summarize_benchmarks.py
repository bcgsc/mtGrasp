#!/usr/bin/env python3

# This script generates a summary of the benchmark outputs (i. runtime in seconds ii. peak memory in MB)
# Input requirements: 1. text file with the path to each benchmark folder 2. csv output file

import sys
import os
import os.path
import pandas as pd

if len(sys.argv[1:]) != 2:
    print("Usage: %s <text filewith the path to each benchmark folder> <output csv file>" % sys.argv[0])
    sys.exit(1)


input = sys.argv[1]
output = sys.argv[2]

file = open(input, 'r')
lines = file.read().splitlines()


step = []
assembly = []
run_time = []
peak_memory = []
library = []


for line in lines:
    directory = "%s/benchmark"%(line)
    
    dir_exists=os.path.isdir(directory)
    if dir_exists == True:
        for file in os.listdir(directory):
            library.append(os.path.basename(line))
            basename = os.path.basename(file)
            assembly.append(basename.split('.')[0])
            step.append(basename.split('.')[1])
            df_bm = pd.read_csv("%s/%s"% (directory, basename), sep="\t")
            run_time.append(float(df_bm.at[0, 's']))
            peak_memory.append(float(df_bm.at[0, 'max_rss']))

df = pd.DataFrame({'Library':library, 'Pipeline Step':step, 'Run time (s)':run_time, 'Peak memory (MB)':peak_memory, 'k/kc':assembly})
df.to_csv(output)

            





