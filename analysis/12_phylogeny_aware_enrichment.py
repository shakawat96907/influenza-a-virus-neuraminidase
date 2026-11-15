#!/usr/bin/env python3
import os, math
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.contingency_tables import StratifiedTable
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo, SeqIO

def compute_or_ci(a, b, c, d):
    # Haldane-Anscombe correction if any cell is 0
    aa, bb, cc, dd = (a + 0.5, b + 0.5, c + 0.5, d + 0.5) if min(a, b, c, d) == 0 else (a, b, c, d)
    odds = (aa * dd) / (bb * cc)
    se = math.sqrt(1/aa + 1/bb + 1/cc + 1/dd)
    ci_low = math.exp(math.log(odds) - 1.96 * se)
    ci_high = math.exp(math.log(odds) + 1.96 * se)
    return odds, ci_low, ci_high

def get_declustered_ids(tree_path, min_dist=0.003):
    tree = Phylo.read(tree_path, "newick")
    kept = list()
    for clade in tree.find_clades(order="level"):
        if clade.is_terminal():
            if not kept:
                kept.append(clade)
            else:
                # Approximate tip-to-tip distance via shared depth
                d_cur = tree.distance(tree.root, clade)
                too_close = False
                for k in kept[-30:]:
                    if abs(d_cur - tree.distance(tree.root, k)) < min_dist and tree.distance(clade, k) < min_dist:
                        too_close = True
                        break
                if not too_close:
                    kept.append(clade)
    return set(t.name for t in kept)

def main():
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, yr_col = col_map["accession"], col_map["group"], col_map["year"]
    meta = meta[meta[acc_col].isin(aln_recs)].copy()
    meta["Is_Asia"] = meta[grp_col].str.strip().str.lower() == "asia"
    meta["Year_Int"] = pd.to_numeric(meta[yr_col], errors="coerce").astype(int)

    asia_ids = set(meta[meta["Is_Asia"]][acc_col])
    nonasia_ids = set(meta[~meta["Is_Asia"]][acc_col])
    n_asia, n_non = len(asia_ids), len(nonasia_ids)

    # Deduplicated protein set (1 per unique protein sequence)
    seen_seqs = set()
    dedup_ids = set()
    for acc in meta[acc_col]:
        s = aln_recs[acc].replace("-", "")
        if s not in seen_seqs:
            seen_seqs.add(s)
            dedup_ids.add(acc)

    # Clade-declustered set (1 representative per close phylogenetic cluster <0.003 subs/site)
    declust_ids = get_declustered_ids("results/trees/tree_nuc_clean_rooted.treefile", min_dist=0.003) & set(meta[acc_col])
    print(f"Dataset Sizes -> Full: {len(meta)} (Asia={n_asia}, Non-Asia={n_non}) | Dedup: {len(dedup_ids)} | Clade-Declustered: {len(declust_ids)}")

    # Load T6 ancestral origins
    t8 = pd.read_csv("results/tables/Supplementary Table 8 - (ancestral_reconstruction_and_homoplasy).csv").set_index("N1_Position")
    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)

    n1_pos = 0
    rows = list()
    years = sorted(meta["Year_Int"].unique())

    for aln_idx, cons_aa in enumerate(cons_gapped):
        if cons_aa == "-":
            continue
        n1_pos += 1
        var_counts = Counter(aln_recs[sid][aln_idx] for sid in aln_recs if aln_recs[sid][aln_idx] not in ("-", cons_aa))
        for var_aa, tot_cnt in var_counts.items():
            glob_pct = (tot_cnt / len(meta)) * 100
            if glob_pct < 1.0:
                continue
            a_has = sum(1 for sid in asia_ids if aln_recs[sid][aln_idx] == var_aa)
            na_has = sum(1 for sid in nonasia_ids if aln_recs[sid][aln_idx] == var_aa)
            a_no, na_no = n_asia - a_has, n_non - na_has

            odds, ci_lo, ci_hi = compute_or_ci(a_has, a_no, na_has, na_no)
            _, p_raw = fisher_exact([[a_has, a_no], [na_has, na_no]], alternative="two-sided")

            # Year-stratified Cochran-Mantel-Haenszel (CMH) test
            tables = list()
            for yr in years:
                sub_a = meta[(meta["Year_Int"] == yr) & (meta["Is_Asia"])][acc_col]
                sub_na = meta[(meta["Year_Int"] == yr) & (~meta["Is_Asia"])][acc_col]
                ya = sum(1 for sid in sub_a if aln_recs[sid][aln_idx] == var_aa)
                yna = sum(1 for sid in sub_na if aln_recs[sid][aln_idx] == var_aa)
                if len(sub_a) > 0 and len(sub_na) > 0:
                    tables.append(np.array([[ya + 0.5, len(sub_a) - ya + 0.5], [yna + 0.5, len(sub_na) - yna + 0.5]]))
            st = StratifiedTable(tables)
            cmh_or = float(st.oddsratio_pooled)
            cmh_lo, cmh_hi = st.oddsratio_pooled_confint()
            cmh_p = float(st.test_null_odds(correction=True).pvalue)

            # Clade-declustered sensitivity test
            dc_a = [sid for sid in declust_ids if sid in asia_ids]
            dc_na = [sid for sid in declust_ids if sid in nonasia_ids]
            dc_ahas = sum(1 for sid in dc_a if aln_recs[sid][aln_idx] == var_aa)
            dc_nahas = sum(1 for sid in dc_na if aln_recs[sid][aln_idx] == var_aa)
            dc_or, _, _ = compute_or_ci(dc_ahas, len(dc_a) - dc_ahas, dc_nahas, len(dc_na) - dc_nahas)
            _, dc_p = fisher_exact([[dc_ahas, len(dc_a) - dc_ahas], [dc_nahas, len(dc_na) - dc_nahas]])

            # Deduplicated sensitivity test
            dd_a = [sid for sid in dedup_ids if sid in asia_ids]
            dd_na = [sid for sid in dedup_ids if sid in nonasia_ids]
            dd_ahas = sum(1 for sid in dd_a if aln_recs[sid][aln_idx] == var_aa)
            dd_nahas = sum(1 for sid in dd_na if aln_recs[sid][aln_idx] == var_aa)
            dd_or, _, _ = compute_or_ci(dd_ahas, len(dd_a) - dd_ahas, dd_nahas, len(dd_na) - dd_nahas)
            _, dd_p = fisher_exact([[dd_ahas, len(dd_a) - dd_ahas], [dd_nahas, len(dd_na) - dd_nahas]])

            ind_orig = int(t8.loc[n1_pos, "Site_Independent_Origins"]) if n1_pos in t8.index else 0
            rows.append({
                "N1_Position": n1_pos, "Variant": f"{cons_aa}{n1_pos}{var_aa}",
                "Global_Pct": round(glob_pct, 2),
                "Asian_Count": a_has, "Asian_Pct": round(a_has / n_asia * 100, 2),
                "NonAsian_Count": na_has, "NonAsian_Pct": round(na_has / n_non * 100, 2),
                "Fisher_OR": round(odds, 3), "Fisher_95CI": f"{ci_lo:.2f}–{ci_hi:.2f}",
                "CI_Low": ci_lo, "CI_High": ci_hi, "Fisher_Raw_P": p_raw,
                "CMH_OR": round(cmh_or, 3), "CMH_95CI": f"{cmh_lo:.2f}–{cmh_hi:.2f}",
                "CMH_CI_Low": cmh_lo, "CMH_CI_High": cmh_hi, "CMH_Raw_P": cmh_p,
                "Declustered_OR": round(dc_or, 3), "Declustered_P": round(dc_p, 4),
                "Dedup_OR": round(dd_or, 3), "Dedup_P": round(dd_p, 4),
                "Independent_Origins": ind_orig
            })

    res = pd.DataFrame(rows)
    res["Fisher_FDR_P"] = multipletests(res["Fisher_Raw_P"], method="fdr_bh")[1]
    res["CMH_FDR_P"] = multipletests(res["CMH_Raw_P"], method="fdr_bh")[1]
    res["Direction"] = np.where(res["Fisher_OR"] > 1, "Enriched in Asia", "Enriched in Non-Asia")
    res = res.sort_values("Fisher_FDR_P")

    res.drop(columns=["CI_Low", "CI_High", "CMH_CI_Low", "CMH_CI_High"]).to_csv(
        "results/tables/Supplementary Table 9 - (regional_enrichment_all_variants).csv", index=False
    )

    sig_main = res[(res["Global_Pct"] >= 3.0) & (res["Fisher_FDR_P"] < 0.05)].copy()
    main_t7 = sig_main[[
        "Variant", "Asian_Pct", "NonAsian_Pct", "Fisher_OR", "Fisher_95CI",
        "Fisher_FDR_P", "CMH_OR", "CMH_95CI", "CMH_FDR_P", "Declustered_P", "Dedup_P", "Independent_Origins"
    ]].copy()
    for c in ["Fisher_FDR_P", "CMH_FDR_P"]:
        main_t7[c] = main_t7[c].apply(lambda x: f"{x:.2e}")
    main_t7.to_csv("results/tables/Main_Table_7_Regional_Enrichment.csv", index=False)

    print("\n=================================================================")
    print(f" MAIN TABLE 7: FDR-Significant Variants (Global >=3%, FDR < 0.05): {len(main_t7)} variants")
    print("=================================================================")
    print(main_t7.to_string(index=False))

    # Render Clean, Overlap-Free Figure 3 (Legends Outside Axes)
    plot_df = res[res["Global_Pct"] >= 3.0].sort_values("Fisher_OR", ascending=True).tail(18)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.5, 7.2), gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.28})

    y_pos = np.arange(len(plot_df))
    log_or = np.log2(plot_df["Fisher_OR"])
    err_lo = log_or - np.log2(plot_df["CI_Low"])
    err_hi = np.log2(plot_df["CI_High"]) - log_or

    cmh_log = np.log2(plot_df["CMH_OR"])
    cmh_elo = cmh_log - np.log2(plot_df["CMH_CI_Low"])
    cmh_ehi = np.log2(plot_df["CMH_CI_High"]) - cmh_log

    ax1.axvline(0, color="#333333", linestyle="--", lw=1.0, zorder=1)
    h1 = ax1.errorbar(log_or, y_pos + 0.14, xerr=[err_lo, err_hi], fmt="o", color="#C0392B", ecolor="#C0392B", elinewidth=1.3, capsize=2.5, markersize=5.5, label="Unstratified Fisher OR (95% CI)")
    h2 = ax1.errorbar(cmh_log, y_pos - 0.14, xerr=[cmh_elo, cmh_ehi], fmt="s", color="#1F618D", ecolor="#1F618D", elinewidth=1.3, capsize=2.5, markersize=5.0, label="Year-Stratified CMH OR (95% CI)")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(plot_df["Variant"], fontsize=9, fontweight="bold")
    ax1.set_xlabel("log₂ Odds Ratio (Non-Asia Enriched  ←  0  →  Asia Enriched)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A   Regional Enrichment Effect Sizes (Unstratified vs. Year-Stratified)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax1.legend(handles=[h1, h2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8, frameon=False)
    ax1.grid(axis="x", linestyle=":", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Asian vs Non-Asian Prevalence Scatter sized by Independent Origins
    all_v = res.copy()
    sig_mask = all_v["Fisher_FDR_P"] < 0.05
    s1 = ax2.scatter(all_v.loc[~sig_mask, "NonAsian_Pct"], all_v.loc[~sig_mask, "Asian_Pct"], s=all_v.loc[~sig_mask, "Independent_Origins"] * 0.45 + 18, color="#95A5A6", alpha=0.65, edgecolors="#555555", lw=0.5, label="Non-significant (FDR ≥ 0.05)")
    s2 = ax2.scatter(all_v.loc[sig_mask, "NonAsian_Pct"], all_v.loc[sig_mask, "Asian_Pct"], s=all_v.loc[sig_mask, "Independent_Origins"] * 0.45 + 22, color="#D62728", alpha=0.85, edgecolors="#111111", lw=0.7, label="FDR-significant (FDR < 0.05)")
    ax2.plot([0, 60], [0, 60], color="#333333", linestyle="--", lw=1.0)

    for idx_i, (_, r) in enumerate(all_v[sig_mask].head(8).iterrows()):
        dy = 2.2 if idx_i % 2 == 0 else -3.2
        ax2.annotate(r["Variant"], xy=(r["NonAsian_Pct"], r["Asian_Pct"]), xytext=(r["NonAsian_Pct"] + 1.5, r["Asian_Pct"] + dy), fontsize=8, fontweight="bold", arrowprops=dict(arrowstyle="-", lw=0.6, color="#444444"))

    ax2.set_xlabel("Non-Asian Control Prevalence (%)", fontsize=9.5, fontweight="bold")
    ax2.set_ylabel("Asian Cohort Prevalence (%)", fontsize=9.5, fontweight="bold")
    ax2.set_title("B   Regional Prevalence & Phylogenetic Recurrence (Size = Indep. Origins)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax2.set_xlim(-2, 60)
    ax2.set_ylim(-2, 60)
    ax2.legend(handles=[s1, s2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8, frameon=False)
    ax2.grid(True, linestyle=":", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.11)
    fig.savefig("results/figures/Figure3_Phylogeny_Aware_Regional_Enrichment_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure3_Phylogeny_Aware_Regional_Enrichment_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved New Figure 3 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
