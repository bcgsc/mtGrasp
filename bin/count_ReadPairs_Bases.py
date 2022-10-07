from Bio import SeqIO
import sys
import os
import shlex
import subprocess
import pandas as pd


library_list=sys.argv[1] 
csv=sys.argv[2] 

file = open(library_list, 'r')
lines = file.read().splitlines()
n_rp = []
n_bases = []
library_id = []
species = []


for line in lines:
        
        directory = "%s"%(line)
        for subdir, dirs, files in os.walk(directory):
            for f in files:
                if os.path.basename(f) == "r1.fq.gz":
                    cmd = 'gunzip -c %s | wc -l'%(os.path.join(subdir,f))
                    
                    run_cmd = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
                    count = int(float(run_cmd.stdout.read()))
                    n_rp.append(int(count/4))
                    n_bases.append(int(count*2*150))
                    library_id.append(line.split('_')[-1])
                    species.append(" ".join(line.split('_')[0:2]))
pd.DataFrame({'Library ID':library_id, 'Species':species, 'Number of read pairs':n_rp, 'Number of bases':n_bases}).to_csv(csv, index=False)
        
