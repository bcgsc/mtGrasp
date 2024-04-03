#!/usr/bin/env python3
'''
mtGrasp: de novo assembly of reference-grade animal mitochondrial genomes
'''
# Wrapper script for mtGrasp
import argparse
import subprocess
import shlex
import sys

MTGRASP_VERSION = 'mtGrasp v1.1.2' # Make sure to edit the version for future releases

parser = argparse.ArgumentParser(description='mtGrasp: de novo assembly of reference-grade animal mitochondrial genomes')
parser.add_argument('-r1', '--read1', help='Forward read fastq.gz file [Required]')
parser.add_argument('-r2', '--read2', help='Reverse read fastq.gz file [Required]')
parser.add_argument('-o', '--out_dir', help='Output directory [Required]')
parser.add_argument('-m', '--mt_gen', help='Mitochondrial genetic code [Required]')
parser.add_argument('-r', '--ref_path', help='Path to the reference fasta file [Required]')

parser.add_argument('-t', '--threads', help='Number of threads [8]', default = 8, type=int)
parser.add_argument('-k', '--kmer', help='k-mer size used in ABySS de novo assembly [91]', default = 91, type=int)
parser.add_argument('-c', '--kc', help='kc [3]', default = 3, type=int)
parser.add_argument('-p', '--gap_filling_p', help='Merge at most N alternate paths during sealer gap filling step [5]',
                    default = 5, type=int)
parser.add_argument('-b', '--sealer_k', help='k-mer size used in sealer gap filling [60,80,100,120]',
                    default = '60,80,100,120')
parser.add_argument('-sf', '--end_recov_sealer_fpr',
                    help='False positive rate for the Bloom filter used by sealer during flanking end recovery [0.01]',
                    default = 0.01, type=float)
parser.add_argument('-sk', '--end_recov_sealer_k',
                    help='k-mer size used in sealer flanking end recovery [60,80,100,120]',
                    default = '60,80,100,120')
parser.add_argument('-i', '--end_recov_p',
                    help='Merge at most N alternate paths during sealer flanking end recovery [5]',
                    default = 5, type=int)
parser.add_argument('-ma', '--mismatch_allowed',
                    help='Maximum number of mismatches allowed in overlaps between the two ends of the mitochondrial assembly [1]',
                    default = 1, type=int)
parser.add_argument('-sub', '--subsample', help='Subsample N read pairs from two paired FASTQ files [2000000]',
                    default = 2000000, type=int)
parser.add_argument('-nsub', '--nosubsample',
                    help='Run mtGrasp using the entire read dataset without subsampling [False]', action='store_true')
parser.add_argument('-a', '--abyss_fpr', help='False positive rate for the Bloom filter used by ABySS [0.005]',
                    default = 0.005, type=float)
parser.add_argument('-s', '--sealer_fpr',
                    help='False positive rate for the Bloom filter used by Sealer during gap filling [0.01]',
                    default = 0.01, type=float)
parser.add_argument('-mp', '--mitos_path', help='Complete path to runmitos.py', default = None)
parser.add_argument('-an', '--annotate', help='Run gene annotation on the final assembly output [False]',
                    action='store_true')
parser.add_argument('-d', '--delete',
                    help='Delete intermediate subdirectories/files once mtGrasp reaches completion [False]',
                    action='store_true')
parser.add_argument('-v', '--version', action="version", version=MTGRASP_VERSION)
parser.add_argument('-u', '--unlock', help='Remove a lock implemented by snakemake on the working directory',
                    action='store_true')
parser.add_argument('-n', '--dry_run', help='Dry-run pipeline', action='store_true')
parser.add_argument('-test', '--test_run', help='Test run mtGrasp to ensure all required dependencies are installed',
                    action='store_true')



args = parser.parse_args()
r1 = args.read1
r2 = args.read2
out_dir = args.out_dir
mt_gen = args.mt_gen
threads = args.threads
dry_run=args.dry_run
kmer = args.kmer
kc = args.kc
ref_path = args.ref_path
sealer_fpr = args.sealer_fpr
abyss_fpr = args.abyss_fpr
p = args.gap_filling_p
sealer_k = args.sealer_k
end_recov_sealer_fpr = args.end_recov_sealer_fpr
end_recov_p = args.end_recov_p
end_recov_sealer_k = args.end_recov_sealer_k
unlock = args.unlock
mismatch_allowed = args.mismatch_allowed
subsample = args.subsample
nosubsample = args.nosubsample
annotate = args.annotate
delete = args.delete
mitos_path = args.mitos_path
test_run = args.test_run

if not test_run and (not r1 or not r2 or not out_dir or not mt_gen or not ref_path):
    parser.error("the following arguments are required: -r1/--read1, -r2/--read2," \
                " -o/--out_dir, -m/--mt_gen, -r/--ref_path")

# If the mitos path isn't supplied, check if 'runmitos.py' is in PATH
if not args.mitos_path:
    output = subprocess.Popen(['which', 'runmitos.py'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, _ = output.communicate()
    path_output = stdout.decode().strip()
    if '/runmitos.py' in path_output:
        mitos_path = path_output.replace('/runmitos.py','')
    else:
        print("Please ensure runmitos.py is available on your PATH or specify the path using -mp")
        sys.exit(1)


# get the directory of the mtgrasp.smk script
string = subprocess.check_output(['which', 'mtgrasp.smk'])
string = string.decode('utf-8')
# remove new line character
string = string.strip()
# split string by '/'
script_dir = '/'.join(string.split('/')[0:-1])

# get the base names of the read files
if not test_run:
    read1_base = r1.split('/')[-1].split('.')[0]
    read2_base = r2.split('/')[-1].split('.')[0]

if test_run:
    subprocess.run(shlex.split(f"{script_dir}/test/test.sh {mitos_path}"), check=True)
elif dry_run:
    subprocess.run(shlex.split(f"snakemake -s {script_dir}/mtgrasp.smk -np --config r1={r1} r2={r2} out_dir={out_dir} "\
                                f"mt_code={mt_gen} k={kmer} kc={kc} ref_path={ref_path} threads={threads} " \
                                f"abyss_fpr={abyss_fpr} sealer_fpr={sealer_fpr} p={p}  sealer_k={sealer_k} " \
                                f"end_recov_sealer_fpr={end_recov_sealer_fpr} end_recov_p={end_recov_p} " \
                                f"end_recov_sealer_k={end_recov_sealer_k} mismatch_allowed={mismatch_allowed} "\
                                f"annotate='N/A' mitos_path={mitos_path} "),
                   check=True)
elif unlock:
    subprocess.run(shlex.split(f"snakemake -s {script_dir}/mtgrasp.smk --unlock --config r1={r1} r2={r2} " \
                                f"out_dir={out_dir} mt_code={mt_gen} k={kmer} kc={kc} ref_path={ref_path} " \
                                f"threads={threads} abyss_fpr={abyss_fpr} sealer_fpr={sealer_fpr} p={p} " \
                                f" sealer_k={sealer_k}  end_recov_sealer_fpr={end_recov_sealer_fpr} " \
                                f"end_recov_p={end_recov_p} end_recov_sealer_k={end_recov_sealer_k} " \
                                f"mismatch_allowed={mismatch_allowed} annotate='N/A' mitos_path={mitos_path} "),
                   check=True)
# if subsample is specified, run the pipeline on the subsampled reads
elif nosubsample and annotate:
    print('Subsampling skipped')
# Run mtgrasp.smk on the entire read dataset without subsampling
    subprocess.run(shlex.split(f"snakemake -s {script_dir}/mtgrasp.smk --cores {threads} -p -k --config r1={r1} " \
                                f"r2={r2} out_dir={out_dir} mt_code={mt_gen} k={kmer} kc={kc} ref_path={ref_path} " \
                                f"threads={threads} abyss_fpr={abyss_fpr} sealer_fpr={sealer_fpr} p={p} " \
                                f" sealer_k={sealer_k}  end_recov_sealer_fpr={end_recov_sealer_fpr} " \
                                f"end_recov_p={end_recov_p} end_recov_sealer_k={end_recov_sealer_k} " \
                                f"mismatch_allowed={mismatch_allowed} annotate='Yes' mitos_path={mitos_path} "),
                   check=True)
elif nosubsample and not annotate:
    print('Subsampling skipped')
# Run mtgrasp.smk on the entire read dataset without subsampling
    subprocess.run(shlex.split(f"snakemake -s {script_dir}/mtgrasp.smk --cores {threads} -p -k --config r1={r1} " \
                                f"r2={r2} out_dir={out_dir} mt_code={mt_gen} k={kmer} kc={kc} ref_path={ref_path} " \
                                f" threads={threads} abyss_fpr={abyss_fpr} sealer_fpr={sealer_fpr} p={p} " \
                                f" sealer_k={sealer_k}  end_recov_sealer_fpr={end_recov_sealer_fpr} " \
                                f"end_recov_p={end_recov_p} end_recov_sealer_k={end_recov_sealer_k} " \
                                f"mismatch_allowed={mismatch_allowed} annotate='No' mitos_path={mitos_path} "),
                   check=True)
elif not nosubsample and annotate :
# By default, run mtgrasp.smk on the subsampled reads
    print('Start subsampling reads')
    subprocess.run(shlex.split(f"sub_then_run_mtgrasp.sh {out_dir} {r1} {r2} {subsample} \
                               {read1_base} {read2_base} {script_dir} {threads} {mt_gen} {kmer} \
                               {kc} {ref_path}  {abyss_fpr} {sealer_fpr} {p} {sealer_k} \
                               {end_recov_sealer_fpr} {end_recov_p} {end_recov_sealer_k} {mismatch_allowed} \
                                'Yes' {mitos_path}"),
                   check=True)
elif not nosubsample and not annotate:
# By default, run mtgrasp.smk on the subsampled reads
    print('Start subsampling reads')
    subprocess.run(shlex.split(f"sub_then_run_mtgrasp.sh {out_dir} {r1} {r2} {subsample} \
                               {read1_base} {read2_base} {script_dir} {threads} {mt_gen} {kmer} \
                               {kc} {ref_path}  {abyss_fpr} {sealer_fpr} {p} {sealer_k} \
                               {end_recov_sealer_fpr} {end_recov_p} {end_recov_sealer_k} {mismatch_allowed} 'No' {mitos_path}"),
                   check=True)

else:
    print('Please double check mtGrasp usage information')
    sys.exit(1)


if delete:
    subprocess.run(shlex.split(f'bash {script_dir}/cleanup.sh {out_dir}'),
                   check=True)
