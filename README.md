# Molecular Evolution and Antigenic Variation of Influenza A(H1N1) Neuraminidase in Asia (2015–2025)

**NCBI Taxonomy ID:** `11320` | **Target Protein:** Neuraminidase (NA, Segment 6)

## Overview

Influenza A virus is a major human respiratory pathogen characterized by continuous genetic and antigenic evolution. Neuraminidase (NA), one of its two major surface glycoproteins, plays an essential role in viral release, motility through respiratory mucus, and host-cell interactions, while serving as the primary target of neuraminidase inhibitor (NAI) antivirals.

This repository provides the reproducible bioinformatics pipeline, sequence curation scripts, curated metadata, codon-aware phylogenetic reconstructions, molecular-clock outlier audits, geographic migration analyses, ancestral state reconstructions, and structural residue numbering audits supporting our comparative evolutionary investigation of human influenza A(H1N1)pdm09 neuraminidase in Asia (2015–2025) alongside a year-matched global control dataset.

## Repository Structure

```text
influenza-a-virus-neuraminidase/
├── CITATION.cff                             # Citation metadata for automated Zenodo DOI archiving
├── run_all.sh                               # End-to-end master reproduction script (Scripts 01–22)
├── README.md                              # Project overview, directory guide, and reproducibility instructions
├── environment.yml                        # Conda environment specification
├── .gitignore                             # Git ignore rules for large/temporary intermediate files
├── data/
├── README.md                          # Documentation of raw and processed sequence datasets
│   ├── raw/                               # Raw sequence records directory
│   └── processed/                         # Curated non-redundant datasets (N = 1,449)
│       ├── final_clean_meta.csv           # Complete metadata & accession list (N = 1,449; Supplementary Table 4)
│       ├── final_clean_prot.fasta         # Curated full-length NA protein sequences (1,449)
│       └── final_clean_nuc.fasta          # Curated frame-matched NA coding nucleotide sequences (1,449)
├── scripts/
├── README.md                          # Execution order and description of data curation scripts
│   └── *.py                               # Retrieval, quality control, deduplication & stratified subsampling
├── analysis/
│   ├── 00_data_curation/                    # Step 0 (T1): NCBI Virus retrieval, QC, deduplication & stratified subsampling
├── README.md                          # Step-by-step documentation of evolutionary & structural analyses
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
│   ├── 11b_clean_fig2_and_tables56.py     # Main Tables 5–6 export & publication Figure 2A–C rendering
│   ├── 12_phylogeny_aware_enrichment.py       # Regional enrichment (Fisher, CMH, clade-declustered; Table 7)
│   ├── 12b_fix_fig3.py                        # 1000 DPI Figure 3A–B rendering with phylogeny-aware overlays
│   ├── 13_shannon_entropy_and_fig4.py         # Shannon entropy & 3-tier variability sensitivity (Figure 4, Table 8)
│   ├── 14_temporal_dynamics_and_fig5.py       # 4-period temporal dynamics & rarefaction (Figure 5, Table 9)
│   ├── 15_glycosylation_and_fig6.py           # Predicted N-linked glycosylation sequon profiling (Figure 6, Table 10)
│   ├── 16_selection_and_fig7.py               # Codon-level HyPhy FEL, MEME & FUBAR selection (Figure 7, Table 11)
│   ├── 17_tetramer_sasa_and_fig8.py           # Biological tetramer (4B7R) SASA & normalized RSA (Main Table 12)
│   ├── 17b_compile_figure8_abc.py             # Single-row composite Figure 8A–C & 4-subunit PyMOL 3D rendering
│   ├── 18_season_matched_vaccine_fig9.py       # Season-matched WHO vaccine divergence (2015–2025; Main Table 13)
│   ├── 18b_render_fig9_fig10_1000dpi.py       # 1000 DPI Figure 9 & PyMOL 3D homotetramer Figure 10 rendering
│   ├── 19_nai_resistance_and_fig11.py         # Canonical NAI resistance screening & H275Y/I223 audit (Figure 11, Table 14)
│   ├── 20_haplotype_cooccurrence_fig12.py     # Multi-variant haplotype greedy set-cover & co-occurrence network (Figure 12, Table 15)
│   ├── 21_epitope_hotspot_scoring_fig13.py    # IEDB/PDB antibody footprint mapping & reproducible hotspot scoring (Figure 13, Table 16)
│   └── 22_weight_sensitivity_catalytic_final.py # Hotspot weight sensitivity & catalytic active-site conservation (Table 17, Supp Fig S2)
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

## Completed Analytical Workflow (`T1`–`T18`)

1. **Dataset Curation & Stratified Subsampling (`analysis/00_data_curation/`):**
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

9. **Temporal Stratification & Equal-Size Rarefaction (`analysis/14_temporal_rarefaction_and_fig5.py`):**
   * Evaluated temporal diversity across 3-period, 4-period (isolating the `2021–2022` COVID-19 bottleneck), and annual (`2015–2025`) schemes using `1,000` random rarefaction iterations (`n = 225`, `n = 36`, `n = 20`) and minimum recurrence filtering (`≥2` sequences), confirming a robust post-2020 evolutionary rebound (**Figure 5A–C**, **Main Table 9**, `Supplementary Table 11`).

10. **Predicted N-Linked Glycosylation Sequon Profiling (`analysis/15_glycosylation_and_fig6.py`, `analysis/15b_clean_fig6_table10.py`):**
   * Scanned canonical `N-X-S/T` (`X ≠ P`) motifs across all 1,449 sequences (`2015–2025`), identifying 7 persistent core sequons (`≥80%` across all periods) and a progressive decline at `N50` (`97.38%` to `78.35%`) driven by `S52N` (**Figure 6A–B**, **Main Table 10**, `Supplementary Table 12`).

11. **Site-Specific Evolutionary Selection Pressure (`analysis/16_selection_and_fig7.py`):**
   * Executed codon-level HyPhy `FEL`, `MEME`, and `FUBAR` analyses across `1,449` frame-matched coding sequences, confirming pervasive purifying selection (`mean dN/dS = 0.258`) and identifying `6` diversifying sites (`p < 0.05` / `PP ≥ 0.90`, including `N200S` and `V264I`) alongside `4` episodic selection sites (**Figure 7A–C**, **Main Table 11**, `Supplementary Table 13`).

12. **Biological Homotetramer SASA/RSA & 3D Structural Mapping (`analysis/17_tetramer_sasa_and_fig8.py`, `analysis/17b_compile_figure8_abc.py`):**
   * Computed monomer vs. biological homotetramer (`PDB 4B7R`) `SASA` and normalized `RSA` (`Tien et al., 2013`), identifying 59 interfacial residues (`ΔSASA ≥ 10 Å²`, including `V453` and `T452`) and compiling single-row Figure 8A–C with PyMOL 3D stick/sphere element-colored views (**Figure 8A–C**, **Main Table 12**, `Supplementary Table 14`).

13. **Season-Matched WHO Vaccine Divergence Profiling (`analysis/18_season_matched_vaccine_fig9.py`, `analysis/18b_render_fig9_fig10_1000dpi.py`):**
    └── 19_nai_resistance_and_fig11.py            # Canonical NAI resistance screening & H275Y/I223 phylogenetic audit (Figure 11, Table 14)
   * Audited `QFQ15359.1` (non-pdm09 historical H1N1) and evaluated percentage divergence from season-matched WHO A(H1N1)pdm09 vaccine virus NA sequences (`2015–2025`), eliminating all `257` legacy `≥50%` artifacts and mapping the `16` top divergent head residues (`≥15%`) onto the 3D homotetramer (**Figures 9 & 10**, **Main Table 13**, `Supplementary Table 15`).

14. **Neuraminidase Inhibitor (`NAI`) Resistance Screening & Phylogenetic Audit (`analysis/19_nai_resistance_and_fig11.py`):**
   * Screened 16 canonical N1 resistance substitutions (`Screened vs. Detected`), disaggregated position `223` (`I223V` detected `n = 8`; `I223R/K/T` = `0`), and audited all `14` `H275Y` isolates across country, year, and global phylogeny (`11` independent phylogenetic origins; **Figure 11A–B**, **Main Table 14**, `Supplementary Table 16`).

15. **Multi-Variant Haplotype Greedy Set-Cover & Co-Occurrence Network (`analysis/20_haplotype_cooccurrence_fig12.py`):**
   * Reconstructed per-isolate haplotypes across Asia-enriched and clade-sweep variants (`n = 750`), showing that `0/750` isolates carry all 7 Asia-enriched mutations simultaneously (`5` strains needed for `≥90%` coverage) due to segregation into two mutually exclusive temporal clades (**Figure 12A–B**, **Main Table 15**, `Supplementary Table 17`).

16. **Experimentally Verified Epitope Footprint Cross-Mapping & Reproducible Hotspot Scoring (`analysis/21_epitope_hotspot_scoring_fig13.py`):**
   * Cross-referenced all `469` canonical N1 residues against curated X-ray crystal structure antibody footprints (`PDB 4QNP`, `6Q23`, `6N6B`, `6PZW`, `7U4E`, `8G3O`; with permanent `DOIs`) and `IEDB` epitopes, separating verified antibody contacts (`n = 41`) from predicted surface candidates (`RSA ≥ 0.25`; `n = 73`) and computing a deterministic `10`-point multi-evidence antigenic hotspot score (**Figure 13A–B**, **Main Table 16**, `Supplementary Tables 18 & 18b`).

17. **Hotspot Weight Sensitivity & Canonical Active-Site Conservation (`analysis/22_weight_sensitivity_catalytic_final.py`):**
   * Verified rank stability between evidence-weighted (`3:2:2:2:1`) and equal-weighted (`2:2:2:2:2`) hotspot scores (`Spearman ρ = 0.9346`, `P < 10⁻¹⁰⁰`; `100%` [`12/12`] top-hotspot overlap; `Supplementary Table 19`, `Supplementary Figure S2`) and confirmed `99.03%–100.00%` conservation across all `20` canonical N1 catalytic triad and framework residues (**Main Table 17**).

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

**Affiliation:** Department of Biochemistry and Molecular Biology, Shahjalal University of Science and Technology (SUST), Sylhet-3114, Bangladesh
