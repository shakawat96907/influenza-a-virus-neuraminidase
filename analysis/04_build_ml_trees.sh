#!/usr/bin/env bash
# Step 04 (Pipeline Step 7 / T3): Maximum-Likelihood Phylogenetics (Nucleotide & Amino-Acid)
set -euo pipefail

mkdir -p results/trees

echo "=========================================================="
echo " [1/2] Inferring Global Nucleotide ML Tree (New Figure 1) "
echo "=========================================================="
iqtree -s results/codon_aln_nuc_cap20.fasta \
       -m MFP \
       -B 1000 \
       -alrt 1000 \
       -T AUTO \
       --seed 42 \
       --prefix results/trees/tree_nuc_cap20

echo "=========================================================="
echo " [2/2] Inferring Global Amino-Acid ML Tree (Supplement)   "
echo "=========================================================="
iqtree -s results/aln_prot_cap20.fasta \
       -m MFP \
       -B 1000 \
       -alrt 1000 \
       -T AUTO \
       --seed 42 \
       --prefix results/trees/tree_prot_cap20

echo "Both ML trees inferred successfully!"
