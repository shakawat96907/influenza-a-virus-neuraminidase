# Data Curation & Stratified Subsampling Pipeline (`T1`)

This directory contains the reproducible Python scripts used to retrieve, filter, deduplicate, and curate the **1,449** human influenza A(H1N1)pdm09 neuraminidase (NA) sequences (`750 Asian` and `699 year-matched non-Asian`) analyzed in this study.

## Execution Order
1. **Metadata Filtering & Deduplication:** Filters NCBI Virus GenBank records (`2015–2025`, human host, `1,380–1,500 nt`, complete CDS, one record per isolate name).
2. **Country-Year Stratified Subsampling (`seed = 42`):** Caps sampling at ≤20 sequences per country per year (`757` Asian and `705` year-matched non-Asian candidates).
3. **CDS Translation & Quality Control:** Verifies intact start codons (`ATG`), absence of internal stop codons, length divisible by 3, `≤1%` ambiguous bases, and protein length `464–474 aa` (`1,456` sequences).
4. **Post-Hoc Consensus Identity & Clock Outlier Filtering (`analysis/07_clean_outliers_and_gard.py`):** Excludes seven non-pdm09/swine-origin variant spillovers or misannotated long-branch outliers (`Supplementary Table 2`), producing the final curated dataset (`data/processed/final_clean_meta.csv`, `N = 1,449`).
