# Molecular Evolution and Antigenic Variation of Influenza A(H1N1) Neuraminidase in Asia (2015–2025)

**NCBI Taxonomy ID:** `11320` | **Target Protein:** Neuraminidase (NA, Segment 6)

## Overview

Influenza A virus is a major human respiratory pathogen characterized by continuous genetic and antigenic evolution. Neuraminidase (NA), one of its two major surface glycoproteins, plays an essential role in viral release, motility through respiratory mucus, and host-cell interactions, while serving as the primary target of neuraminidase inhibitor (NAI) antivirals.

This repository provides the reproducible bioinformatics pipeline, sequence curation scripts, curated metadata, codon-aware phylogenetic reconstructions, molecular-clock outlier audits, geographic migration analyses, ancestral state reconstructions, and structural residue numbering audits supporting our comparative evolutionary investigation of human influenza A(H1N1)pdm09 neuraminidase in Asia (2015–2025) alongside a year-matched global control dataset.

## Repository Structure

```text
influenza-a-virus-neuraminidase/
├── README.md                              # Project overview, directory guide, and reproducibility instructions
├── environment.yml                        # Conda environment specification
├── run_all.sh                             # Master script to execute the completed analysis pipeline
├── .gitignore                             # Git ignore rules for large/temporary intermediate files
├── data/
│   ├── README.md                          # Documentation of raw and processed sequence datasets
│   ├── raw/                               # Raw sequence records directory
│   └── processed/                         # Curated non-redundant datasets (N = 1,449)
│       ├── final_clean_meta.csv           # Complete metadata & accession list (N = 1,449; Supplementary Table 4)
│       ├── final_clean_prot.fasta         # Curated full-length NA protein sequences (1,449)
│       └── final_clean_nuc.fasta          # Curated frame-matched NA coding nucleotide sequences (1,449)
├── scripts/
│   ├── README.md                          # Execution order and description of data curation scripts
│   └── *.py                               # Retrieval, quality control, deduplication & stratified subsampling
├── analysis/
│   ├── README.md                          # Step-by-step documentation of evolutionary & structural analyses
│   ├── 01_mafft_and_consensus.sh          # Multiple sequence alignment (MAFFT v7.526) & consensus generation
│   ├── 02_numbering_audit_4B7R.py         # Pairwise structural alignment (PDB 4B7R) & N1/N2 coordinate audit
│   ├── 03_codon_alignment.py              # CDS reading-frame trimming & pal2nal.pl codon alignment
│   ├── 04_build_ml_trees.sh               # IQ-TREE Maximum-Likelihood inference (Nucleotide & Amino-Acid)
│   ├── 05_tree_audit_and_rooting.py       # Model/support comparison & outgroup rooting (KU322051)
│   ├── 06_clock_and_recombination.py      # TreeTime root-to-tip clock regression & branch-length audit
│   ├── 07_clean_outliers_and_gard.py      # Outlier exclusion (N = 1,449), re-alignment & HyPhy GARD
│   ├── 08_export_manuscript_draft.py      # Consensus regeneration, 4B7R audit & manuscript summary tables
│   ├── 09_finalize_trees_and_fig1.py      # Clean tree topology comparison, clock fit & Figure 1A–B rendering
│   ├── 10_treetime_mugration.py           # Discrete geographic mugration & Asia-enriched clade detection
│   ├── 11_ancestral_and_fig2.py           # TreeTime ancestral sequence reconstruction & homoplasy audit
│   └── 11b_clean_fig2_and_tables56.py     # Main Tables 5–6 export & publication Figure 2A–C rendering
└── results/
    ├── aln_prot_clean.fasta               # Clean MAFFT protein alignment (1,449 sequences; 493 columns)
    ├── codon_aln_nuc_clean.fasta          # Clean codon-aware nucleotide alignment (1,449 sequences; 1,479 nt)
    ├── consensus_prot_clean.fasta         # Gapped majority-rule consensus sequence (493 alignment columns)
    ├── consensus_prot_cap20_nogap.fasta     # Canonical ungapped consensus sequence (469 amino acids)
    ├── trees/                             # Rooted & unrooted IQ-TREE Newick files and model logs
    ├── treetime_clock_clean/              # Clean TreeTime root-to-tip clock regression outputs (β = 3.25e-3)
    ├── mugration/                         # TreeTime discrete geographic mugration models & annotated trees
    ├── ancestral/                         # TreeTime ancestral codon/protein reconstructions & annotated tree
    ├── figures/                           # High-resolution publication figures (Figure 1A–B, Fig 2A–C, Fig S1)
    └── tables/                            # Main Tables 5–6, Supplementary Tables 1–8 & outlier audit logs
```

## Completed Analytical Workflow (`T1`–`T9`)

1. **Dataset Curation & Stratified Subsampling (`scripts/`):**
   * Retrieved human A(H1N1)pdm09 NA coding sequences (2015–2025) from NCBI Virus, filtered for complete coding sequences (1,380–1,500 nt), ≤1% ambiguous bases, intact reading frames, and deduplicated isolate names.
   * Applied country- and year-capped stratified sampling (≤20 sequences/country/year) and excluded seven non-pdm09/variant spillover or misannotated long-branch outliers (`Supplementary Table 2`), yielding **750 Asian** and **699 non-Asian** sequences (**total *N* = 1,449**).
2. **Multiple Sequence Alignment & Consensus Construction (`analysis/01_mafft_and_consensus.sh`, `analysis/08_export_manuscript_draft.py`):**
   * Aligned all **1,449** non-redundant NA protein sequences using `MAFFT` (v7.526), generating a 493-column alignment and a canonical **469-amino-acid** ungapped consensus sequence (`H1N1_NA_consensus`).
3. **Structural Alignment & Residue Numbering Audit (`analysis/02_numbering_audit_4B7R.py`, `analysis/08_export_manuscript_draft.py`):**
   * Aligned the 469-aa consensus NA sequence to Chain A of **PDB ID: 4B7R** (positions 83–469), confirming **96.64% sequence identity** (374/387 identical residues) and **100.00% head-domain coverage**, and exported a full N1/PDB/N2 mapping table (`Supplementary Table 3`).
4. **Codon-Aware Global Phylogenetics & Molecular-Clock Outlier Resolution (`analysis/03_codon_alignment.py` – `analysis/09_finalize_trees_and_fig1.py`):**
   * Inferred maximum-likelihood phylogenies in `IQ-TREE` (v3.1.4) with `1,000 UFBoot` and `1,000 SH-aLRT` replicates for both the global nucleotide alignment (`GTR+F+I+R3`; `N = 1,449`; **Figure 1A**) and global amino-acid alignment (`HIVW+I+R3`; `871` unique taxa; **Figure 1B**), rooted on `KU322051` (2015, Italy).
   * Verified molecular-clock coherence in `TreeTime` (v0.12.1), confirming linear substitution accumulation (`β = 3.25 × 10⁻³` substitutions/site/year; inferred root date `2012.0` [`2011.97`]; **Supplementary Figure S1**).
5. **Discrete Geographic Ancestral State Reconstruction (`analysis/10_treetime_mugration.py`):**
   * Reconstructed intercontinental dispersal across six global regions, identifying **338 inter-regional transitions** (**234 exports from Asia** and **49 introductions into Asia**; **Main Table 5**, `Supplementary Table 6`) alongside **66 Asia-enriched monophyletic subclusters** (`≥80%` Asian representation; `Supplementary Table 7`).
6. **Ancestral Sequence Reconstruction & Homoplasy Analysis (`analysis/11_ancestral_and_fig2.py`, `analysis/11b_clean_fig2_and_tables56.py`):**
   * Performed `TreeTime ancestral` reconstruction across the 1,449-taxon phylogeny to quantify phylogenetically independent `parent → child` substitution events across all 469 N1 sites, distinguishing clade-sweeping variants from recurrent homoplastic substitutions (**Figure 2A–C**, **Main Table 6**, `Supplementary Table 8`).

7. **Phylogeny-Aware & Year-Stratified Regional Enrichment (`analysis/12_phylogeny_aware_enrichment.py`, `analysis/12b_fix_fig3.py`):**
   * Evaluated regional enrichment across all 59 variants (≥1.0% global prevalence) using two-sided Fisher's exact tests (`BH FDR`), year-stratified Cochran–Mantel–Haenszel (`CMH`) tests, protein-deduplicated sensitivity analysis (`n = 718`), and clade-declustered phylogeny-pruned sensitivity analysis (`n = 538`) paired with `TreeTime` independent origins (**Figure 3A–B**, **Main Table 7**, `Supplementary Table 9`).

8. **Site-Wise Shannon Entropy & Cut-Off Sensitivity (`analysis/13_shannon_entropy_and_fig4.py`):**
   * Computed information-theoretic Shannon entropy (`H`, bits) and `1 − consensus frequency` across all 469 N1 positions, confirming near-perfect rank concordance (`Spearman ρ = 0.9991`) and nested stability across lenient, baseline, and stringent variability thresholds (**Figure 4A–C**, **Main Table 8**, `Supplementary Table 10`).

## Ongoing Revision Analyses (`T9`–`T18`)

* Phylogeny-aware regional enrichment (`Fisher's exact test + BH FDR` & independent origin comparison) and within-year Shannon entropy / temporal rarefaction (`T7`–`T8`)
* Site-specific selection analysis (`HyPhy MEME / FEL / FUBAR`) and co-occurring mutation / epistasis network (`T9`–`T10`)
* Predicted N-linked glycosylation sequon profiling (`N-X-S/T`) and biological tetramer (`PDB 4B7R`) SASA/RSA structural analysis (`T11`–`T12`)
* Season-matched WHO vaccine divergence (`2015–2025`), antiviral resistance screening, reproducible antigenic scoring, and Zenodo DOI archival release (`T13`–`T18`)

## Reproducibility

To reproduce the completed pipeline locally:

```bash
conda env create -f environment.yml
conda activate na
bash run_all.sh
```

## Authors & Correspondence

* **Md. Shakawat Hossain** — First Author
* **Rubaiya Tasnim** — Second Author
* **Dr. Farjana Islam** — Corresponding Author (`farjana-bmb@sust.edu`)

* **Md. Shakawat Hossain** — First Author
* **Rubaiya Tasnim** — Second Author
* **Dr. Farjana Islam** — Corresponding Author (`farjana-bmb@sust.edu`)
