#!/usr/bin/env python3
'''
#This script incorporates a portion of the "circle_check.py" script from the MitoZ project (Meng et al., 2019) - 
# "circle_check.py" script Link: 
#  https://github.com/linzhi2013/MitoZ/blob/a43e62862dd7d7253ed4ce79b4a3e972ebb75218/version_2.3/useful_scripts/circle_check.py
#This script first checks if the assembly sequence is above 12kbp and has an overlapping region between the 2 flanking regions
#If an overlap is found, "abyss-mergepairs" will be used to merge the 2 end segments containing the overlap region into a single sequence
#In case of non-overlapping regions, two end segments (each comprising 200 base pairs) will be obtained from the initial assembly sequence. 
# These segments will then be utilized to fabricate a "fake gap" structure, such as end(200bp) + NNNNNNNNNN + start(200bp).
'''

import sys
import shlex
import subprocess
import os
from Bio import SeqIO

if len(sys.argv[1:]) < 9:
    print(f"Usage: {sys.argv[0]} <assembly> <bf> <r1> <r2> <out_dir> <threads> <p> <mismatch_allowed> <k> [k ..]")
    sys.exit(1)

assembly = sys.argv[1]
bf = sys.argv[2]
r1 = sys.argv[3]
r2 = sys.argv[4]
out_dir = sys.argv[5]
threads = sys.argv[6]
p = sys.argv[7]
mismatch_allowed = int(sys.argv[8])
k = sys.argv[9:]
k = " ".join(k)


cmd = f'mkdir -p {out_dir}'
cmd_shlex = shlex.split(cmd)
subprocess.call(cmd_shlex)


# Check if the assembly sequence is above 12kbp and has an overlapping region between the 2 flanking regions

def seq_distance(seq1=None, seq2=None):
    "Find the number of mismatches between the sequences"
    if len(seq1) != len(seq2):
        print(seq1)
        print(seq2)
        sys.exit("Error: seq1 and seq2 have different lengths!")

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    mismatch = 0
    for i, base in enumerate(seq1):
        if base != seq2[i]:
            mismatch += 1
    return mismatch

def reverse_complement(sequence):
    "Reverse complement the sequence"
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return ''.join([complement[base] for base in sequence[::-1]])

def create_fake_fastq_read(sequence, end):
    "Create the PE fastq read for end recovery"
    quality_score = "F" # Quality score  for all bases
    if end == 'end':
        name1 = "@sequence1/1"
        with open(f"{out_dir}/end_sequence.fastq", "w", encoding="utf-8") as file1:
        # Write the first sequence to the first file
            file1.write(f"{name1}\n{sequence}\n+\n{quality_score * len(sequence)}\n")
    elif end == 'start':
        name2 = "@sequence2/2"
        with open(f"{out_dir}/start_sequence.fastq", "w", encoding="utf-8") as file2:
        # Write the second sequence to the second file
            reverse_com = reverse_complement(sequence)
            file2.write(f"{name2}\n{reverse_com}\n+\n{quality_score * len(sequence)}\n")
    else:
        sys.exit("Error: end must be defined as either 'end' or 'start'!")



rec_list = []
for rec in SeqIO.parse(assembly, 'fasta'):
    rec_list.append(str(rec.seq))
seq = rec_list[0]
seqlen = len(seq)
if seqlen < 12000:
    print("Assembly sequence is below 12kbp. No need to perform overlap check")
    fh_out = open(f'{out_dir}/fake_gap_unfilled.fa', 'w', encoding="utf-8")
    end_bp = []
    start_bp = []
    end_bp.append(f"{seq[-200:]}")
    start_bp.append(f"{seq[:200]}")
    list_seq = end_bp + ["NNNNNNNNNN"] + start_bp
    string = ''.join(list_seq)
    outline = ">" + rec.id + ' fake_gap'+ "\n" + string
    print(outline, file=fh_out)
    fh_out.close()

potential_seq = seq[12000:]
overlap_len = 5
match_start_pos = 0
tmp = 0
while True:
    subseq_5 = seq[0:overlap_len]
    match_count = potential_seq.count(subseq_5)
    if match_count > 1:
        overlap_len += 1
    elif match_count == 1:
        match_start_pos = 12000 + potential_seq.find(subseq_5)
        break
    else:
        tmp = 1
        break

if tmp == 1:
    print("Overlap not found")
    fh_out = open(f'{out_dir}/fake_gap_unfilled.fa', 'w', encoding="utf-8")
    end_bp = []
    start_bp = []
    end_bp.append(f"{seq[-200:]}")
    start_bp.append(f"{seq[:200]}")
    list_seq = end_bp + ["NNNNNNNNNN"] + start_bp
    string = ''.join(list_seq)
    outline = ">" + rec.id + ' fake_gap'+ "\n" + string
    print(outline, file=fh_out)
    fh_out.close()


start_pos = match_start_pos
stop_pos = match_start_pos + overlap_len
subseq_3 = seq[start_pos:stop_pos]
if seqlen >= 12000 and start_pos !=0 and seq_distance(subseq_5, subseq_3) < mismatch_allowed or seq_distance(subseq_5, subseq_3) == mismatch_allowed:
    print(f"Overlap found ( {start_pos + 1} bp to {stop_pos + 1} bp )")

    end_seq = seq[start_pos:]
    start_seq = seq[0:len(seq[start_pos:])]

    create_fake_fastq_read(end_seq, 'end')
    create_fake_fastq_read(start_seq, 'start')
    cmd = f'mtgrasp_abyss_mergepairs.sh {out_dir} end_sequence.fastq start_sequence.fastq {mismatch_allowed}'
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)
    # check if {outdir}/out_merged.fastq is not empty
    if os.stat(f"{out_dir}/out_merged.fastq").st_size != 0:
        print("abyss_mergepairs merged the ends")
    else:
        fh_out = open(f'{out_dir}/fake_gap_unfilled.fa', 'w', encoding="utf-8")
        print("abyss_mergepairs failed to merge the ends")
        print("Creating fake gap")
        end_bp = []
        start_bp = []
        end_bp.append(f"{seq[-200:]}")
        start_bp.append(f"{seq[:200]}")
        list_seq = end_bp + ["NNNNNNNNNN"] + start_bp
        string = ''.join(list_seq)
        outline = ">" + rec.id + ' fake_gap'+ "\n" + string
        print(outline, file=fh_out)
        fh_out.close()
elif seqlen >= 12000 and seq_distance(subseq_5, subseq_3) > mismatch_allowed:
    fh_out = open(f'{out_dir}/fake_gap_unfilled.fa', 'w', encoding="utf-8")
    print("Overlap exceeds maximum mismatch allowed")
    end_bp = []
    start_bp = []
    end_bp.append(f"{seq[-200:]}")
    start_bp.append(f"{seq[:200]}")
    list_seq = end_bp + ["NNNNNNNNNN"] + start_bp
    string = ''.join(list_seq)
    outline = ">" + rec.id + ' fake_gap'+ "\n" + string
    print(outline, file=fh_out)
    fh_out.close()
else:
    pass


# check if {outdir}/fake_gap_unfilled.fa exists
if os.path.exists(f"{out_dir}/fake_gap_unfilled.fa") and os.stat(f"{out_dir}/fake_gap_unfilled.fa").st_size != 0:
    print("fake_gap_unfilled.fa created")
    print("Start sealer gap filling")
    cmd = f'abyss-sealer -b{bf} -j {threads} -vv {k} -P {p} -o {out_dir}/fake_gap_filled -S {out_dir}/fake_gap_unfilled.fa {r1} {r2}'
    print(cmd)
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)
    print("Start adding flanks back to the original assembly")
    cmd = f'mtgrasp_check_filled_add_flanks.py {out_dir}/fake_gap_filled_log.txt {out_dir}/fake_gap_filled_scaffold.fa {assembly} {out_dir}/flank_added_assembly.fa {out_dir}/fake_gap_unfilled.fa'
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)


elif os.path.exists(f"{out_dir}/out_merged.fastq") and os.stat(f"{out_dir}/out_merged.fastq").st_size != 0:
    print('Start adding abyss-mergepairs merged sequence back to the original assembly')
    for record in SeqIO.parse(f'{out_dir}/out_merged.fastq', "fastq"):
        merged_seq = str(record.seq)
        flanks_extracted_seq = seq[len(seq[start_pos:]): start_pos]
        scenario_1_seq = merged_seq + flanks_extracted_seq
        scenario_2_seq = flanks_extracted_seq + merged_seq

        with open(f'{out_dir}/flank_added_assembly.fa', 'w', encoding="utf-8") as f:
            f.write(">" + "Scenario 1 Sequence" + "\n" + scenario_1_seq + "\n" + ">" + "Scenario 2 Sequence" + "\n" + scenario_2_seq)

else:
    # print error message
    print("Error: No fake_gap_unfilled.fa or out_merged.fastq created")
    