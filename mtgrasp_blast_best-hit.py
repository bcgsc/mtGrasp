#!/usr/bin/env python3
import sys
import shlex
import subprocess
from collections import namedtuple

'''
BLASTs query against the database specified, and parses the output.
Outputs best hit per query (highest coverage + sequence identity)
blastn must be in your PATH
'''

# Helper function -- decide if current hit is better than a previous hit
def is_better_hit(coverage, perc_id, prev_hit):
    if coverage > prev_hit.coverage:
        return True
    if coverage == prev_hit.coverage and perc_id > prev_hit.perc_id:
        return True
    return False

if len(sys.argv[1:]) != 2:
    print("Usage: %s <query> <blast db name>" % sys.argv[0])
    sys.exit(1)

query = sys.argv[1]
blast_dbname = sys.argv[2]

Alignment = namedtuple("Alignment", ['query', 'ref', 'perc_id', 'coverage', 'start_ref', 'end_ref', 'ref_title'])

hits = {} # query -> Alignment

cmd = 'blastn -query %s -num_threads 4 -db %s  -outfmt "6 qaccver saccver length nident qlen qstart qend slen sstart send stitle"' %\
      (query, blast_dbname)
cmd_shlex = shlex.split(cmd)

run_blast = subprocess.Popen(cmd_shlex, stdout=subprocess.PIPE, universal_newlines=True)
for line in iter(run_blast.stdout.readline, ''):
    (query, ref, length_align, nident, qlen, qstart, qend, slen, sstart, send, stitle) = line.strip().split('\t')

    # Calculate percent ID
    nident, qstart, qend = float(nident), int(qstart), int(qend)
    perc_id = nident/(qend - qstart)*100

    # Calculate coverage
    qstart, qend, qlen = int(qstart), int(qend), float(qlen)
    coverage = (qend - qstart)/qlen*100

    if (query in hits and is_better_hit(coverage, perc_id, hits[query])) or \
        query not in hits:
        hits[query] = Alignment(query=query, ref=ref, perc_id=perc_id, coverage=coverage, start_ref=int(sstart), end_ref=int(send), ref_title=stitle)

# Print out the best hits
print("query", "ref", "percent_identity", 'coverage', 'start_ref', 'end_ref', 'ref_title', sep="\t")

for hit in sorted(hits.values(), key=lambda x: (x.ref, x.start_ref)):
    print(hit.query, hit.ref, hit.perc_id, hit.coverage, hit.start_ref, hit.end_ref, hit.ref_title, sep="\t")

