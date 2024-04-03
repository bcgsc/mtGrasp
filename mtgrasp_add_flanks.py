#!/usr/bin/env python3

"""
Script for adding the gap-filled joint flanking sequence to the original assembly which has been trimmed of the flanking regions.
Usage: 
python3 script.py <fasta file with gap-filled construct> <original assembly> <flank size> <output fasta prefix>
"""


from Bio import SeqIO
import sys 

gap_filled_construct = sys.argv[1]
assembly = sys.argv[2]
flank_size = sys.argv[3]
output = sys.argv[4]

# Check for the correct number of command-line arguments
if len(sys.argv[1:]) != 4:
    print("Usage: %s <fasta file that stores the 'fake' gap-filled postSealer DNA construct> <original assembly> <flank size used for creating fake-gap construct> <output fasta file prefix>" % sys.argv[0])
    sys.exit(1)

for construct in SeqIO.parse(gap_filled_construct, "fasta"):
    for record in SeqIO.parse(assembly, "fasta"):
        flanks_extracted_seq = str(record.seq)[int(flank_size):-int(flank_size)] # extract the flanks from the original assembly
        
         # Having two scenarios to make sure the tRNA-phe isn't divided between the two ends, allowing us to have a full tRNA-phe annotation for start site standardization
        scenario_1_seq = str(construct.seq) + flanks_extracted_seq # scenario 1: the gap filled construct is inserted to the start of the original assembly
        scenario_2_seq = flanks_extracted_seq + str(construct.seq) # scenario 2: the gap filled construct is inserted to the end of the original assembly

        # write the two scenarios to the output file
        with open(output, 'w') as f:
               f.write(">" + "Scenario 1 Sequence" + "\n" + scenario_1_seq + "\n" + ">" + "Scenario 2 Sequence" + "\n" + scenario_2_seq)