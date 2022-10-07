import sys 
from Bio import SeqIO

fasta=sys.argv[1]
flank_length=sys.argv[2]
output=sys.argv[3]
end_bp = []
start_bp = []

for record in SeqIO.parse(fasta, "fasta"):
    end_bp.append(">%s\n%s" % (record.description, record.seq[-int(flank_length):]))
    start_bp.append("%s" % (record.seq[:int(flank_length)]))

list = end_bp + ["NNNNNNNNNN"] + start_bp
string = ''.join(list)

fakegap_file = open(output, 'w+')
fakegap_file.write(string)
fakegap_file.close()


