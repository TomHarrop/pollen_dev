#!/usr/bin/env python3

from Bio import SeqIO

fasta_file = snakemake.input['fasta']
seqlen_file = snakemake.output['seqlen']

# read fasta
seqlens = [[x.id, str(len(x))] for x in SeqIO.parse(fasta_file, 'fasta')]

# write seqlength
with open(seqlen_file, 'wt') as f:
    for x in seqlens:
        f.write('\t'.join(x))
        f.write('\n')
