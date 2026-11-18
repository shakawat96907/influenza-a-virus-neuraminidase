# Core Bioinformatics & Evolutionary Analysis Pipeline

1. `01_mafft_and_consensus.sh` — Multiple sequence alignment (MAFFT v7.526 `--auto`) and generation of gapped (493-column) and canonical ungapped (469-aa) consensus sequences.
2. `02_numbering_audit_4B7R.py` — Global pairwise alignment against crystal structure PDB ID: `4B7R` (Chain A; 96.64% sequence identity) and cross-numbering standardization (MAFFT alignment, canonical N1, PDB `4B7R`, and classical N2 numbering).

17. `17_tetramer_sasa_and_fig8.py` — Biological homotetramer (`PDB ID: 4B7R`) vs. monomer `SASA` and size-normalized `RSA` (`Tien et al., 2013`) calculation (**Main Table 12**, `Supplementary Table 14`).
18. `17b_compile_figure8_abc.py` — Single-row composite **Figure 8A–C** rendering (`SASA`, `RSA`, and 4-subunit `CPK` PyMOL 3D homotetramer structure).
19. `18_season_matched_vaccine_fig9.py` — Season-matched WHO A(H1N1)pdm09 vaccine divergence profiling (`2015–2025`) and `QFQ15359.1` provenance audit (**Main Table 13**, `Supplementary Table 15`).
20. `18b_render_fig9_fig10_1000dpi.py` — `1000 DPI` **Figure 9** (`A–B`) and PyMOL 3D homotetramer **Figure 10** (`CPK` element-colored vaccine-divergent landmarks) rendering.
21. `19_nai_resistance_and_fig11.py` — Canonical N1 neuraminidase inhibitor (`NAI`) resistance screening (`Screened vs. Detected`), residue `223` disaggregation, and global phylogenetic clustering audit (**Figure 11A–B**, **Main Table 14**, `Supplementary Table 16`).
