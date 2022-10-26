#!/usr/bin/env python3

from Bio import SeqIO
import sys 

gap_filled_construct = sys.argv[1]
assembly = sys.argv[2]
flank_size = sys.argv[3]
output = sys.argv[4]

if len(sys.argv[1:]) != 4:
    print("Usage: %s <fasta file that stores the 'fake' gap-filled postSealer DNA construct> <original assembly> <flank size used for creating fake-gap construct> <output fasta file prefix>" % sys.argv[0])
    sys.exit(1)

for construct in SeqIO.parse(gap_filled_construct, "fasta"):
    for record in SeqIO.parse(assembly, "fasta"):
        flanks_extracted_seq = str(record.seq)[int(flank_size):-int(flank_size)]
        scenario_1_seq = str(construct.seq) + flanks_extracted_seq
        scenario_2_seq = flanks_extracted_seq + str(construct.seq)

        with open(output, 'w') as f:
               f.write(">" + "Scenario 1 Sequence" + "\n" + scenario_1_seq + "\n" + ">" + "Scenario 2 Sequence" + "\n" + scenario_2_seq)