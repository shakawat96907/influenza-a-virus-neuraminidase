#!/usr/bin/env bash
set -euo pipefail
mkdir -p results

mafft --auto --thread 4 data/processed/final_cap20_prot.fasta > results/aln_prot_cap20.fasta
seqkit stats results/aln_prot_cap20.fasta

python -c "
from collections import Counter
seqs = []
with open('results/aln_prot_cap20.fasta') as f:
    curr = []
    for line in f:
        if line.startswith('>'):
            if curr: seqs.append(''.join(curr))
            curr = []
        else:
            curr.append(line.strip())
    if curr: seqs.append(''.join(curr))

cons = ''.join(Counter(col).most_common(1) for col in zip(*seqs))
with open('results/consensus_prot_cap20.fasta', 'w') as out:
    out.write('>H1N1_NA_consensus\n' + cons + '\n')
"

seqkit seq -g results/consensus_prot_cap20.fasta > results/consensus_prot_cap20_nogap.fasta
seqkit stats results/consensus_prot_cap20_nogap.fasta
