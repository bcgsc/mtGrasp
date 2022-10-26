#!/usr/bin/env python3

# Wrapper script for mtGasp
import argparse
import subprocess
import shlex

parser = argparse.ArgumentParser(description='Usage of mtGasp')
parser.add_argument('-r1', '--read1', help='Forward read fastq.gz file', required=True)
parser.add_argument('-r2', '--read2', help='Reverse read fastq.gz file', required=True)
parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
parser.add_argument('-m', '--mt_gen', help='Mitochondrial genetic code', required=True)
parser.add_argument('-t', '--threads', help='Number of threads', nargs='?', const=1, default = 1)
parser.add_argument('-k', '--kmer', help='k-mer size', nargs='?', const=1, default = 96)
parser.add_argument('-c', '--kc', help='kc', nargs='?', const=1, default = 3)
parser.add_argument('-p', '--ref_path', help='k-mer size', required=True)
parser.add_argument('-n', '--dry_run', help='Dry-run pipeline', action='store_true')


args = parser.parse_args()
r1 = args.read1
r2 = args.read2
out_dir = args.out_dir
mt_gen = args.mt_gen
threads = args.threads
kmer = args.kmer
kc = args.kc
ref_path = args.ref_path
dry_run = args.dry_run 


# get the directory of the mtgasp.smk script
string = subprocess.check_output(['which', 'mtgasp.smk'])
string = string.decode('utf-8')
# remove new line character
string = string.strip()
# split string by '/'
script_dir = '/'.join(string.split('/')[0:-1])

if dry_run:
    subprocess.run(shlex.split(f'snakemake -s {script_dir}/mtgasp.smk -n --config r1={r1} r2={r2} out_dir={out_dir} mt_code={mt_gen} k={kmer} kc={kc} ref_path={ref_path}'))
else:
# Run mtgasp.smk
    subprocess.run(shlex.split(f'snakemake -s {script_dir}/mtgasp.smk --cores {threads} --config r1={r1} r2={r2} out_dir={out_dir} mt_code={mt_gen} k={kmer} kc={kc} ref_path={ref_path}'))





