import sys 
from Bio import SeqIO
import re
 

fasta=sys.argv[1]
list = []

for record in SeqIO.parse(fasta, "fasta"):
  string = str(record.seq)
  result = re.findall(r'ATG[A-Z]{9}CAT', string)
  list.append(result[0])

print(list)

