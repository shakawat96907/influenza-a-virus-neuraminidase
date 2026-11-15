# Molecular Evolution and Antigenic Variation of Influenza A(H1N1) Neuraminidase in Asia (2015–2025)

**NCBI Taxonomy ID:** `11320` | **Target Protein:** Neuraminidase (NA, Segment 6)

## Overview

Influenza A virus is a major human respiratory pathogen characterized by continuous genetic and antigenic evolution. Neuraminidase (NA), one of its two major surface glycoproteins, plays an essential role in viral release, motility through respiratory mucus, and host-cell interactions, while serving as the primary target of neuraminidase inhibitor (NAI) antivirals.

This repository provides the reproducible bioinformatics pipeline, sequence curation scripts, curated metadata, codon-aware phylogenetic reconstructions, molecular-clock outlier audits, geographic migration analyses, and structural residue numbering audits supporting our comparative evolutionary investigation of human influenza A(H1N1)pdm09 neuraminidase in Asia (2015–2025) alongside a year-matched global control dataset.

## Repository Structure

```text
influenza-a-virus-neuraminidase/
├── README.md                          # Project overview, directory guide, and reproducibility instructions
├── environment.yml                    # Conda environment specification
├── run_all.sh                         # Master script to execute the completed analysis pipeline
├── .gitignore                         # Git ignore rules for large/temporary intermediate files
├── data/
│   ├── README.md                      # Documentation of raw and processed sequence datasets
│   ├── raw/                           # Raw sequence records directory
│   └── processed/                     # Curated non-redundant datasets (N = 1,449)
│       ├── final_clean_meta.csv       # Complete metadata & accession list (N = 1,449; Supplementary Table 4)
│       ├── final_clean_prot.fasta     # Curated full-length NA protein sequences (1,449)
│       └── final_clean_nuc.fasta      # Curated frame-matched NA coding nucleotide sequences (1,449)
├── scripts/
│   ├── README.md                      # Execution order and description of data curation scripts
│   └── *.py                           # Retrieval, quality control, deduplication & stratified subsampling
├── analysis/
│   ├── README.md                      # Step-by-step documentation of evolutionary & structural analyses
│   ├── 01_mafft_and_consensus.sh      # Multiple sequence alignment (MAFFT v7.526) & consensus generation
│   ├── 02_numbering_audit_4B7R.py     # Pairwise structural alignment (PDB 4B7R) & N1/N2 coordinate audit
│   ├── 03_codon_alignment.py          # CDS reading-frame trimming & pal2nal.pl codon alignment
│   ├── 04_build_ml_trees.sh           # IQ-TREE Maximum-Likelihood inference (Nucleotide & Amino-Acid)
│   ├── 05_tree_audit_and_rooting.py   # Model/support comparison & outgroup rooting (KU322051)
│   ├── 06_clock_and_recombination.py  # TreeTime root-to-tip clock regression & branch-length audit
│   ├── 07_clean_outliers_and_gard.py  # Outlier exclusion (N = 1,449), re-alignment & HyPhy GARD
│   ├── 08_export_manuscript_draft.py  # Consensus regeneration, 4B7R audit & manuscript summary tables
│   ├── 09_finalize_trees_and_fig1.py  # Clean tree topology comparison, clock fit & Figure 1A–B rendering
│   └── 10_treetime_mugration.py       # Discrete geographic mugration & Asia-enriched clade detection
└── results/
    ├── aln_prot_clean.fasta           # Clean MAFFT protein alignment (1,449 sequences; 493 columns)
    ├── codon_aln_nuc_clean.fasta      # Clean codon-aware nucleotide alignment (1,449 sequences; 1,479 nt)
    ├── consensus_prot_cap20.fasta     # Gapped majority-rule consensus sequence (493 alignment columns)
    ├── consensus_prot_cap20_nogap.fasta # Canonical ungapped consensus sequence (469 amino acids)
    ├── trees/                         # Rooted & unrooted IQ-TREE Newick files and model logs
    ├── treetime_clock_clean/          # Clean TreeTime root-to-tip clock regression outputs (β = 3.25e-3)
    ├── mugration/                     # TreeTime discrete geographic mugration models & annotated trees
    ├── figures/                       # High-resolution publication figures (Figure 1A–B, Fig S1)
    └── tables/                        # Supplementary Tables 1–7 & outlier exclusion audit logs
```

## Completed Analytical Workflow (`T1`–`T5`)

1. **Dataset Curation & Stratified Subsampling (`scripts/`):**
   * Retrieved human A(H1N1)pdm09 NA coding sequences (2015–2025) from NCBI Virus.
   * Filtered for complete coding sequences (1,380–1,500 nt), ≤1% ambiguous bases, and absence of premature stop codons or frameshifts, followed by isolate-level deduplication.
   * Applied a country- and year-capped stratified sampling scheme (maximum 20 sequences per country per year), followed by sequence-identity and phylogenetic outlier auditing (excluding seven divergent non-pdm09/variant spillover or misannotated outliers; `Supplementary Table 2`), yielding a final curated dataset of **750 Asian** and **699 non-Asian** sequences (**total *N* = 1,449**).
2. **Multiple Sequence Alignment & Consensus Construction (`analysis/01_mafft_and_consensus.sh`, `analysis/08_export_manuscript_draft.py`):**
   * Aligned all **1,449** non-redundant NA protein sequences using `MAFFT` (v7.526) (`--auto`, FFT-NS-2 progressive alignment strategy), generating a 493-column alignment and a canonical **469-amino-acid** ungapped consensus sequence (`H1N1_NA_consensus`).
3. **Structural Alignment & Residue Numbering Audit (`analysis/02_numbering_audit_4B7R.py`, `analysis/08_export_manuscript_draft.py`):**
   * Performed global pairwise alignment between the 469-aa consensus NA sequence and Chain A of the reference crystal structure **PDB ID: 4B7R** (387 resolved residues, positions 83–469), confirming **96.64% sequence identity** (374/387 identical residues) and **100.00% coverage** across the crystallographic head domain (82.52% full-length coverage).
   * Standardized all coordinates to canonical **N1 numbering**, validated conserved active-site and antiviral resistance benchmark residues (`E119`, `Q136`, `I223`, `S247`, `H275`, `R293`, `N295`), and exported a complete cross-numbering mapping table (`results/tables/Supplementary Table 3 - (numbering_audit_table).csv`).
4. **Codon-Aware Global Phylogenetics & Molecular-Clock Outlier Resolution (`analysis/03_codon_alignment.py` – `analysis/09_finalize_trees_and_fig1.py`):**
   * Constructed a codon-aware nucleotide alignment (`1,479` nt across `493` codons) via `pal2nal.pl` (v14) and inferred maximum-likelihood phylogenies using `IQ-TREE` (v3.1.4) with `1,000 UFBoot` and `1,000 SH-aLRT` replicates for both the global nucleotide dataset (`GTR+F+I+R3`; `N = 1,449`; **Figure 1A**) and the global amino-acid dataset (`HIVW+I+R3`; `871` unique taxa; **Figure 1B**).
   * Rooted phylogenies on the earliest 2015 reference isolate (`KU322051`, Italy) and verified molecular-clock coherence using `TreeTime` (v0.12.1), confirming linear substitution accumulation (`β = 3.25 × 10⁻³` substitutions/site/year; inferred root date `2012.0` [`2011.97`]; **Supplementary Figure S1**).
5. **Discrete Geographic Ancestral State Reconstruction (`analysis/10_treetime_mugration.py`):**
   * Reconstructed intercontinental migration pathways across six geographic regions (`Asia`, `North America`, `Europe`, `Africa`, `Oceania`, `South America`), identifying **338 inter-regional transitions** (**234 exports from Asia** and **49 introductions into Asia**; `Supplementary Table 6`) alongside **66 Asia-enriched monophyletic subclusters** (`≥80%` Asian representation; `Supplementary Table 7`).

## Ongoing Revision Analyses (`T6`–`T18`)

* Ancestral state reconstruction & phylogenetic trunk mutation tracking (`T6`)
* Phylogeny-aware regional enrichment (`Fisher's exact test + BH FDR`) and within-year Shannon entropy / temporal rarefaction (`T7`–`T8`)
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
