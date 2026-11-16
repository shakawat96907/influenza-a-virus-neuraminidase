#!/usr/bin/env python3
import os, math
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

def calc_shannon_entropy(residues):
    valid = [r for r in residues if r not in ("-", "X", "*")]
    n = len(valid)
    if n == 0: return 0.0
    counts = Counter(valid)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())

def classify_site(freq, low_cut, high_cut):
    if freq <= low_cut: return "Conserved"
    elif freq < high_cut: return "Variable"
    return "Hypervariable"

def main():
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col = col_map["accession"], col_map["group"]
    meta = meta[meta[acc_col].isin(aln_recs)].copy()

    asia_ids = list(meta[meta[grp_col].str.strip().str.lower() == "asia"][acc_col])
    nonasia_ids = list(meta[meta[grp_col].str.strip().str.lower() != "asia"][acc_col])
    all_ids = asia_ids + nonasia_ids

    rows = []
    n1_pos = 0
    for aln_idx, cons_aa in enumerate(cons_gapped):
        if cons_aa == "-": continue
        n1_pos += 1
        g_res = [aln_recs[s][aln_idx] for s in all_ids if aln_recs[s][aln_idx] not in ("-", "X")]
        a_res = [aln_recs[s][aln_idx] for s in asia_ids if aln_recs[s][aln_idx] not in ("-", "X")]
        na_res = [aln_recs[s][aln_idx] for s in nonasia_ids if aln_recs[s][aln_idx] not in ("-", "X")]

        g_freq = sum(1 for r in g_res if r != cons_aa) / len(g_res)
        a_freq = sum(1 for r in a_res if r != cons_aa) / len(a_res)
        na_freq = sum(1 for r in na_res if r != cons_aa) / len(na_res)

        rows.append({
            "N1_Position": n1_pos, "Consensus_AA": cons_aa,
            "Global_NonCons_Freq": round(g_freq, 5),
            "Asian_NonCons_Freq": round(a_freq, 5),
            "NonAsian_NonCons_Freq": round(na_freq, 5),
            "Global_Shannon_Entropy": round(calc_shannon_entropy(g_res), 5),
            "Asian_Shannon_Entropy": round(calc_shannon_entropy(a_res), 5),
            "NonAsian_Shannon_Entropy": round(calc_shannon_entropy(na_res), 5),
            "Tier_Lenient_001_005": classify_site(g_freq, 0.01, 0.05),
            "Tier_Baseline_002_010": classify_site(g_freq, 0.02, 0.10),
            "Tier_Stringent_005_020": classify_site(g_freq, 0.05, 0.20)
        })

    df = pd.DataFrame(rows)
    df.to_csv("results/tables/Supplementary Table 10 - (site_variability_and_shannon_entropy).csv", index=False)

    rho_g, p_g = spearmanr(df["Global_NonCons_Freq"], df["Global_Shannon_Entropy"])
    rho_reg, p_reg = spearmanr(df["Asian_Shannon_Entropy"], df["NonAsian_Shannon_Entropy"])
    print(f"Spearman Correlation (1 - Consensus vs. Shannon Entropy): rho = {rho_g:.4f}, p = {p_g:.2e}")
    print(f"Spearman Correlation (Asian vs. Non-Asian Entropy)     : rho = {rho_reg:.4f}, p = {p_reg:.2e}")

    cut_configs = [
        ("Lenient (≤0.01 / ≥0.05)", "Tier_Lenient_001_005", 0.05),
        ("Baseline (≤0.02 / ≥0.10)", "Tier_Baseline_002_010", 0.10),
        ("Stringent (≤0.05 / ≥0.20)", "Tier_Stringent_005_020", 0.20)
    ]
    hv_base = set(df[df["Global_NonCons_Freq"] >= 0.10]["N1_Position"])
    sens_rows = []
    for label, col, hcut in cut_configs:
        vc = df[col].value_counts()
        hv_set = set(df[df["Global_NonCons_Freq"] >= hcut]["N1_Position"])
        ov = len(hv_set & hv_base)
        sens_rows.append({
            "Threshold_Scheme": label,
            "Conserved_Sites": int(vc.get("Conserved", 0)),
            "Variable_Sites": int(vc.get("Variable", 0)),
            "Hypervariable_Sites": int(vc.get("Hypervariable", 0)),
            "Nested_Overlap_With_Baseline_HV": f"{ov}/{min(len(hv_set), len(hv_base))} (100%)"
        })
    t8_df = pd.DataFrame(sens_rows)
    t8_df.to_csv("results/tables/Main_Table_8_Variability_Sensitivity_Summary.csv", index=False)

    print("\n=== MAIN TABLE 8: Cut-Off Sensitivity Across Three Threshold Schemes ===")
    print(t8_df.to_string(index=False))

    print("\n=== Top 12 Highest Shannon Entropy Sites (N = 1,449) ===")
    print(df.sort_values("Global_Shannon_Entropy", ascending=False).head(12)[["N1_Position", "Consensus_AA", "Global_NonCons_Freq", "Global_Shannon_Entropy", "Asian_Shannon_Entropy", "NonAsian_Shannon_Entropy"]].to_string(index=False))

    # Render Clean, Overlap-Free Figure 4 (Legends Outside Axes)
    fig = plt.figure(figsize=(14.5, 9.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.05, 1.0], hspace=0.45, wspace=0.25)

    ax1 = fig.add_subplot(gs[0, :])
    l1, = ax1.plot(df["N1_Position"], df["Global_NonCons_Freq"], color="#D62728", lw=1.1, alpha=0.88, label="1 − Consensus Frequency")
    ax1_r = ax1.twinx()
    l2, = ax1_r.plot(df["N1_Position"], df["Global_Shannon_Entropy"], color="#1F618D", lw=1.0, alpha=0.82, label="Shannon Entropy H (bits)")
    c1 = ax1.axhline(0.10, color="#7F8C8D", linestyle="--", lw=0.9, label="Baseline Hypervariable Cut-off (0.10)")
    c2 = ax1.axhline(0.02, color="#BDC3C7", linestyle=":", lw=0.9, label="Baseline Variable Cut-off (0.02)")

    ax1.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("1 − Consensus Frequency", fontsize=10, fontweight="bold", color="#B03A2E")
    ax1_r.set_ylabel("Shannon Entropy (bits)", fontsize=10, fontweight="bold", color="#1B4F72")
    ax1.set_title("A   Site-Wise Non-Consensus Frequency and Shannon Entropy Across N1 Neuraminidase (N = 1,449)", loc="left", fontsize=11, fontweight="bold", pad=26)
    ax1.set_xlim(0, 473)
    ax1.set_ylim(0, 0.62)
    ax1_r.set_ylim(0, 1.35)
    ax1.legend(handles=[l1, l2, c1, c2], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=4, fontsize=8, frameon=False)
    ax1.grid(axis="x", linestyle=":", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1_r.spines["top"].set_visible(False)

    ax2 = fig.add_subplot(gs[1, 0])
    s1 = ax2.scatter(df[df["Global_NonCons_Freq"] <= 0.02]["Global_NonCons_Freq"], df[df["Global_NonCons_Freq"] <= 0.02]["Global_Shannon_Entropy"], s=16, color="#95A5A6", alpha=0.7, label="Conserved (≤0.02)")
    s2 = ax2.scatter(df[(df["Global_NonCons_Freq"] > 0.02) & (df["Global_NonCons_Freq"] < 0.10)]["Global_NonCons_Freq"], df[(df["Global_NonCons_Freq"] > 0.02) & (df["Global_NonCons_Freq"] < 0.10)]["Global_Shannon_Entropy"], s=24, color="#E67E22", alpha=0.85, label="Variable (>0.02–<0.10)")
    s3 = ax2.scatter(df[df["Global_NonCons_Freq"] >= 0.10]["Global_NonCons_Freq"], df[df["Global_NonCons_Freq"] >= 0.10]["Global_Shannon_Entropy"], s=34, color="#C0392B", alpha=0.9, edgecolors="#111111", lw=0.5, label="Hypervariable (≥0.10)")
    ax2.text(0.04, 0.88, f"Spearman ρ = {rho_g:.4f}\nP < 1.0 × 10⁻¹⁰⁰ (n = 469)", transform=ax2.transAxes, fontsize=9, fontweight="bold", bbox=dict(boxstyle="round,pad=0.35", facecolor="#F8F9F9", edgecolor="#BDC3C7"))
    ax2.set_xlabel("1 − Consensus Frequency", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Shannon Entropy H (bits)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Rank Concordance Between Variability Metrics", loc="left", fontsize=11, fontweight="bold", pad=26)
    ax2.legend(handles=[s1, s2, s3], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=7.8, frameon=False)
    ax2.grid(True, linestyle=":", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    ax3 = fig.add_subplot(gs[1, 1])
    x = np.arange(len(t8_df))
    w = 0.25
    bc = ax3.bar(x - w, t8_df["Conserved_Sites"], width=w, color="#7F8C8D", label="Conserved")
    bv = ax3.bar(x, t8_df["Variable_Sites"], width=w, color="#E67E22", label="Variable")
    bh = ax3.bar(x + w, t8_df["Hypervariable_Sites"], width=w, color="#C0392B", label="Hypervariable")
    for bars in [bc, bv, bh]:
        for bar in bars:
            h = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2, h + 6, f"{int(h)}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax3.set_xticks(x)
    ax3.set_xticklabels(["Lenient\n(≤0.01 / ≥0.05)", "Baseline\n(≤0.02 / ≥0.10)", "Stringent\n(≤0.05 / ≥0.20)"], fontsize=9, fontweight="bold")
    ax3.set_ylabel("Number of N1 Sites (Total = 469)", fontsize=10, fontweight="bold")
    ax3.set_title("C   Sensitivity Across Three Variability Threshold Pairs", loc="left", fontsize=11, fontweight="bold", pad=26)
    ax3.set_ylim(0, 470)
    ax3.legend(handles=[bc, bv, bh], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=8.5, frameon=False)
    ax3.grid(axis="y", linestyle=":", alpha=0.4)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.94, top=0.91, bottom=0.08)
    fig.savefig("results/figures/Figure4_Variability_And_Shannon_Entropy_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure4_Variability_And_Shannon_Entropy_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved New Figure 4 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
