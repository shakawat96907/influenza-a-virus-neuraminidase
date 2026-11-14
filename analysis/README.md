# Core Bioinformatics & Evolutionary Analysis Pipeline

1. `01_mafft_and_consensus.sh` — Multiple sequence alignment (MAFFT v7.526 `--auto`) and generation of gapped (493-column) and canonical ungapped (469-aa) consensus sequences.
2. `02_numbering_audit_4B7R.py` — Global pairwise alignment against crystal structure PDB ID: `4B7R` (Chain A; 96.64% sequence identity) and cross-numbering standardization (MAFFT alignment, canonical N1, PDB `4B7R`, and classical N2 numbering).
