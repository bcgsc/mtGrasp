#!/usr/bin/env python3

import sys
from Bio import SeqIO
usage = """
Description
    This script is a modified version of the circle_check.py script from the MitoZ project (Meng et al., 2019)
    Link: https://github.com/linzhi2013/MitoZ/blob/a43e62862dd7d7253ed4ce79b4a3e972ebb75218/version_2.3/useful_scripts/circle_check.py
    Approach:
    - Checking whether the sequences are circular when the sequences have
    - length >= 12Kbp
    - Create a fake gap construct based on the length of overlapping region 
Usage
    python3 {0}  <in.fasta> <mismatch_allowed> <output>
""".format(sys.argv[0])

if len(sys.argv) != 4:
    print(usage)
    sys.exit(0)

in_f, mismatch_allowed,output = sys.argv[1:4]
mismatch_allowed = int(mismatch_allowed)

def seq_distance(seq1=None, seq2=None):
    if len(seq1) != len(seq2):
        print(seq1)
        print(seq2)
        sys.exit("Error: seq1 and seq2 have different lengths!")

    seq1 = seq1.upper()
    seq2 = seq2.upper()

    mismatch = 0
    for i in range(0, len(seq1)):
        if seq1[i] != seq2[i]:
            mismatch += 1
    return (mismatch)

fh_out = open(output, 'w')


for rec in SeqIO.parse(in_f, 'fasta'):
    seq = str(rec.seq)
    seqlen = len(seq)
    if seqlen < 12000:
        end_bp = []
        start_bp = []
        end_bp.append("%s" % (seq[-200:]))
        start_bp.append("%s" % (seq[:200]))
        list = end_bp + ["NNNNNNNNNN"] + start_bp
        string = ''.join(list)
        outline = ">" + rec.id + ' fake_gap'+ "\n" + string
        print(outline, file=fh_out)
        continue
    potential_seq = seq[12000:]
    overlap_len = 5
    match_start_pos = 0
    tmp = 0
    while(1):
        subseq_5 = seq[0:overlap_len]
        match_count = potential_seq.count(subseq_5)
        if match_count > 1:
            overlap_len += 1
        elif match_count == 1:
            match_start_pos = 12000 + potential_seq.find(subseq_5)
            break
        else:
            print("not found!")
            tmp = 1
            break

    if tmp == 1:
        print("A: not found!")
        end_bp = []
        start_bp = []
        end_bp.append("%s" % (seq[-200:]))
        start_bp.append("%s" % (seq[:200]))
        list = end_bp + ["NNNNNNNNNN"] + start_bp
        string = ''.join(list)
        outline = ">" + rec.id + ' fake_gap'+ "\n" + string
        print(outline, file=fh_out)
        continue

    start_pos = match_start_pos
    stop_pos = match_start_pos + overlap_len
    subseq_3 = seq[start_pos:stop_pos]
    if seq_distance(subseq_5, subseq_3) < mismatch_allowed:
        end_bp = []
        start_bp = []

        end_bp.append("%s" % (seq[start_pos:]))
        start_bp.append("%s" % (seq[0:len(seq[start_pos:])]))

        list = end_bp + ["NNNNNNNNNN"] + start_bp
        string = ''.join(list)
        outline = ">" + rec.id + ' fake_gap'+ "\n" + string
        print(outline, file=fh_out)   
        break
    else:
        print("A: not found!")
        end_bp = []
        start_bp = []
        end_bp.append("%s" % (seq[-200:]))
        start_bp.append("%s" % (seq[:200]))
        list = end_bp + ["NNNNNNNNNN"] + start_bp
        string = ''.join(list)
        outline = ">" + rec.id + ' fake_gap'+ "\n" + string
        print(outline, file=fh_out)
        break 
        
            

fh_out.close()
