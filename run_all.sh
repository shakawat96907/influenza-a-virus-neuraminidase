#!/usr/bin/env bash
set -euo pipefail
echo "=== Running Complete A(H1N1)pdm09 NA Revision Pipeline (Scripts 01–22) ==="
bash analysis/01_mafft_and_consensus.sh
python3 analysis/02_numbering_audit_4B7R.py
python3 analysis/03_codon_alignment.py
bash analysis/04_build_ml_trees.sh
python3 analysis/05_tree_audit_and_rooting.py
python3 analysis/06_clock_and_recombination.py
python3 analysis/07_clean_outliers_and_gard.py
python3 analysis/09_finalize_trees_and_fig1.py
python3 analysis/10_treetime_mugration.py
python3 analysis/11_ancestral_and_fig2.py
python3 analysis/11b_clean_fig2_and_tables56.py
python3 analysis/12_phylogeny_aware_enrichment.py
python3 analysis/12b_fix_fig3.py
python3 analysis/13_shannon_entropy_and_fig4.py
python3 analysis/14_temporal_dynamics_and_fig5.py
python3 analysis/15_glycosylation_and_fig6.py
python3 analysis/16_selection_and_fig7.py
python3 analysis/17_tetramer_sasa_and_fig8.py
python3 analysis/17b_compile_figure8_abc.py
python3 analysis/18_season_matched_vaccine_fig9.py
python3 analysis/18b_render_fig9_fig10_1000dpi.py
python3 analysis/19_nai_resistance_and_fig11.py
python3 analysis/20_haplotype_cooccurrence_fig12.py
python3 analysis/21_epitope_hotspot_scoring_fig13.py
python3 analysis/22_weight_sensitivity_catalytic_final.py
echo "=== All 22 analytical scripts completed successfully! ==="
