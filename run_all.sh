#!/usr/bin/env bash
set -euo pipefail

# Master execution script for completed analyses (T1 - T5; N = 1,449)
bash analysis/01_mafft_and_consensus.sh
python analysis/02_numbering_audit_4B7R.py
python analysis/03_codon_alignment.py
bash analysis/04_build_ml_trees.sh
python analysis/05_tree_audit_and_rooting.py
python analysis/06_clock_and_recombination.py
python analysis/07_clean_outliers_and_gard.py
python analysis/08_export_manuscript_draft.py
python analysis/09_finalize_trees_and_fig1.py
python analysis/10_treetime_mugration.py
