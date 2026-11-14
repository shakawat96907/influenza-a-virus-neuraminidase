# Molecular Evolution and Antigenic Variation of Influenza A(H1N1) Neuraminidase in Asia (2015–2025)

**NCBI Taxonomy ID:** 11320 | **Target Protein:** Neuraminidase (NA)

## Overview
Influenza A virus is a major human respiratory pathogen characterized by continuous genetic and antigenic evolution. Neuraminidase (NA), one of its two major surface glycoproteins, plays an essential role in viral release, motility through respiratory mucus, and host-cell interactions, while serving as the primary target of neuraminidase inhibitor (NAI) antivirals.

This repository provides the reproducible bioinformatics pipeline, sequence curation scripts, curated metadata, multiple sequence alignment workflows, and structural residue numbering audits supporting our comparative evolutionary investigation of human influenza A(H1N1)pdm09 neuraminidase in Asia (2015–2025) alongside a year-matched global control dataset.

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
│   └── processed/                         # Curated non-redundant datasets (N = 1,456)
│       ├── final_cap20_meta.csv           # Complete metadata & accession list (Supplementary Table S1)
│       ├── final_cap20_prot.fasta         # Curated full-length NA protein sequences (1,456)
│       ├── final_cap20_nuc.fasta          # Curated full-length NA coding nucleotide sequences (1,456)
│       └── final_cap20_year_table.csv     # Temporal distribution summary across 2015–2025
├── scripts/
│   ├── README.md                          # Execution order and description of data curation scripts
│   └── *.py                               # Retrieval, quality control, deduplication & stratified subsampling scripts
├── analysis/
│   ├── README.md                          # Step-by-step documentation of evolutionary & structural analyses
│   ├── 01_mafft_and_consensus.sh          # Multiple sequence alignment (MAFFT v7.526) & consensus generation
│   └── 02_numbering_audit_4B7R.py         # Pairwise structural alignment (PDB 4B7R) & N1/N2 coordinate audit
└── results/
    ├── consensus_prot_cap20.fasta         # Gapped majority-rule consensus sequence (493 alignment columns)
    ├── consensus_prot_cap20_nogap.fasta   # Canonical ungapped consensus sequence (469 amino acids)
    ├── figures/                           # Publication figures directory
    └── tables/
        └── numbering_audit_table.csv      # Cross-numbering conversion table (Supplementary Table S3)
```

## Completed Analytical Workflow
1. **Dataset Curation & Stratified Subsampling (`scripts/`):**
   - Retrieved human A(H1N1)pdm09 NA coding sequences (2015–2025) from NCBI Virus.
   - Filtered for complete coding sequences (1,380–1,500 nt), ≤1% ambiguous bases, and absence of premature stop codons or frameshifts, followed by isolate-level deduplication.
   - Applied a country- and year-capped stratified sampling scheme (maximum 20 sequences per country per year), yielding a balanced dataset of **754 Asian** and **702 non-Asian** sequences (**total *N* = 1,456**).
2. **Multiple Sequence Alignment & Consensus Construction (`analysis/01_mafft_and_consensus.sh`):**
   - Aligned all 1,456 non-redundant NA protein sequences using **MAFFT (v7.526)** (`--auto`, FFT-NS-2 progressive alignment strategy), generating a 493-column alignment and a canonical **469-amino-acid** ungapped consensus sequence (`H1N1_NA_consensus`).
3. **Structural Alignment & Residue Numbering Audit (`analysis/02_numbering_audit_4B7R.py`):**
   - Performed global pairwise alignment between the 469-aa consensus NA sequence and Chain A of the reference crystal structure **PDB ID: 4B7R** (387 resolved residues, positions 83–469), confirming **96.64% sequence identity** (374/387 identical residues) and **100.00% coverage** across the crystallographic head domain (**82.52%** full-length coverage).
   - Standardized all coordinates to canonical **N1 numbering**, validated conserved active-site and antiviral resistance benchmark residues (`E119`, `Q136`, `I223`, `S247`, `H275`, `R293`, `N295`), and exported a complete cross-numbering mapping table (`results/tables/numbering_audit_table.csv`).

## Ongoing Revision Analyses
- Codon-aware nucleotide alignment and global maximum-likelihood phylogeny (IQ-TREE)
- Root-to-tip clock outlier detection (TreeTime) and recombination screening (HyPhy GARD)
- Geographic migration reconstruction (TreeTime mugration) and ancestral state reconstruction
- Phylogeny-aware regional enrichment, Shannon entropy, temporal rarefaction, N-linked glycosylation sequon profiling, biological tetramer SASA/RSA, season-matched WHO vaccine divergence, antiviral resistance screening, and reproducible antigenic scoring.

## Reproducibility
To reproduce the completed pipeline locally:
```bash
conda env create -f environment.yml
conda activate na
bash run_all.sh
```

## Authors & Correspondence
- **Md. Shakawat Hossain** — First Author
- **Rubaiya Tasnim** — Second Author
- **Dr. Farjana Islam** — Corresponding Author (`farjana-bmb@sust.edu`)

Department of Biochemistry and Molecular Biology  
Shahjalal University of Science and Technology (SUST)  
Sylhet-3114, Bangladesh
