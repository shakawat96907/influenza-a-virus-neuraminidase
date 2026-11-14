# Data Retrieval, Quality Control, and Stratified Subsampling Scripts

This directory contains the 15 Python scripts used for NCBI Virus sequence retrieval, quality control filtering, isolate-level deduplication, and country-year capped stratified sampling:

1. `download_cds.py` & `download_extra_cds.py` — Retrieve complete H1N1 NA coding sequences (CDS) and metadata.
2. `qc_and_filter.py` & `qc_extra.py` — Filter for human host, complete CDS length (1,380–1,500 nt), ≤1% ambiguous bases, and absence of internal stop codons/frameshifts.
3. `populate_meta.py` & `filter_and_select.py` — Parse metadata, assign geographic regions, and deduplicate by isolate name.
4. `get_missing_accessions.py` & `find_missing_accessions.py` — Identify and retrieve missing sequence accessions.
5. `combine_nonasia.py` & `analyze_nonasia.py` — Consolidate and summarize non-Asian global control sequences.
6. `control_matching.py` & `control_matching_human_cap20.py` — Perform year-matched, country-capped stratified subsampling (max 20 sequences/country/year).
7. `rebuild_cap20.py` & `build_final_cap20.py` — Generate final balanced FASTA and metadata files (`N = 1,456`).
8. `verify_analysis.py` — Validate sequence counts and metadata consistency across all strata.
