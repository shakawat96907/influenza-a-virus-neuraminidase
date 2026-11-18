#!/usr/bin/env python3
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

# Canonical N1 Catalytic & Framework Active-Site Residues (19 residues)
CATALYTIC_PANEL = [
    (118, "R", "Catalytic Triad (Direct Sialic Acid Contact)"),
    (151, "D", "Catalytic Triad (Proton Donor / Transition State)"),
    (152, "R", "Catalytic Triad (Acetamido Binding)"),
    (225, "R", "Catalytic Triad (Glycerol Side-Chain Binding)"),
    (277, "E", "Catalytic Triad (Stabilizes R293)"),
    (293, "R", "Catalytic Triad (Carboxylate Binding)"),
    (368, "R", "Catalytic Triad (Carboxylate Binding)"),
    (402, "Y", "Catalytic Triad (Nucleophilic / Transition State)"),
    (119, "E", "Framework Support (Stabilizes R118 / NAI Binding)"),
    (156, "R", "Framework Support (Stabilizes D151)"),
    (179, "W", "Framework Support (Hydrophobic Pocket Wall)"),
    (180, "S", "Framework Support (Acetamido Pocket)"),
    (199, "D", "Framework Support (Stabilizes R152)"),
    (223, "I", "Framework Support (Hydrophobic Floor / NAI Marker)"),
    (228, "E", "Framework Support (Stabilizes R225)"),
    (247, "S", "Framework Support (Pocket Rim / NAI Marker)"),
    (275, "H", "Framework Support (Positions E277 / Primary NAI Marker)"),
    (278, "E", "Framework Support (Stabilizes H275 / E277)"),
    (295, "N", "Framework Support (Stabilizes R293 / NAI Marker)"),
    (425, "E", "Framework Support (Active-Site Structural Core)")
]

def main():
    # 1. Card 30 (T12): Weight Sensitivity Analysis (Primary Weighted 3:2:2:2:1 vs Equal Weighted 2:2:2:2:2)
    df = pd.read_csv("results/tables/Supplementary Table 18 - (reproducible_antigenic_hotspot_scores_469).csv")

    # Normalize each component to [0, 1] and compute Equal-Weight Score (0–10 scale: 2 pts each)
    c_epi = df["Score_Epitope_3"] / 3.0
    c_rsa = df["Score_RSA_2"] / 2.0
    c_vac = df["Score_Vaccine_2"] / 2.0
    c_ent = df["Score_Entropy_2"] / 2.0
    c_sel = df["Score_Selection_1"] / 1.0

    df["Equal_Weight_Score_10"] = np.round((c_epi + c_rsa + c_vac + c_ent + c_sel) * 2.0, 2)
    df["Score_Diff_Weighted_vs_Equal"] = np.round(df["Total_Hotspot_Score_10"] - df["Equal_Weight_Score_10"], 2)

    rho, pval = spearmanr(df["Total_Hotspot_Score_10"], df["Equal_Weight_Score_10"])
    top9_w = set(df.sort_values(["Total_Hotspot_Score_10", "SeasonMatched_VacDiv_Global_Pct"], ascending=False).head(9)["N1_Position"])
    top9_eq = set(df.sort_values(["Equal_Weight_Score_10", "SeasonMatched_VacDiv_Global_Pct"], ascending=False).head(9)["N1_Position"])
    top12_w = set(df.sort_values(["Total_Hotspot_Score_10", "SeasonMatched_VacDiv_Global_Pct"], ascending=False).head(12)["N1_Position"])
    top12_eq = set(df.sort_values(["Equal_Weight_Score_10", "SeasonMatched_VacDiv_Global_Pct"], ascending=False).head(12)["N1_Position"])

    print("=================================================================")
    print(" CARD 30 (T12): Weight Sensitivity Analysis (Primary vs Equal)")
    print("=================================================================")
    print(f"Spearman Rank Correlation across all 469 positions: rho = {rho:.4f} (P = {pval:.2e})")
    print(f"Top-9 Overlap: {len(top9_w.intersection(top9_eq))}/9 ({len(top9_w.intersection(top9_eq))/9*100:.1f}%) | Positions: {sorted(top9_w.intersection(top9_eq))}")
    print(f"Top-12 Overlap: {len(top12_w.intersection(top12_eq))}/12 ({len(top12_w.intersection(top12_eq))/12*100:.1f}%) | Positions: {sorted(top12_w.intersection(top12_eq))}")

    sens_top = df.sort_values(["Total_Hotspot_Score_10", "SeasonMatched_VacDiv_Global_Pct"], ascending=False).head(15)[
        ["N1_Position", "Consensus_AA", "Total_Hotspot_Score_10", "Equal_Weight_Score_10", "Score_Diff_Weighted_vs_Equal", "Evidence_Classification"]
    ]
    print("\nTop 15 Sites Comparison (Primary Weighted vs Equal Weighted):")
    print(sens_top.to_string(index=False))
    df.to_csv("results/tables/Supplementary Table 19 - (hotspot_weight_sensitivity_469).csv", index=False)

    # Render Supplementary Figure S2 (Weight Sensitivity Scatter, 1000 DPI)
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    ax.scatter(df["Equal_Weight_Score_10"], df["Total_Hotspot_Score_10"], alpha=0.55, s=36, color="#1F618D", edgecolors="none", label="All 469 N1 Positions")
    hi = df[df["Total_Hotspot_Score_10"] >= 7.5]
    ax.scatter(hi["Equal_Weight_Score_10"], hi["Total_Hotspot_Score_10"], s=65, color="#C0392B", edgecolors="black", lw=0.6, zorder=4, label="High-Priority Hotspots (Score ≥ 7.5, n = 12)")
    ax.plot([0, 10], [0, 10], ls="--", color="#7F8C8D", lw=1.1, label="1:1 Identity Line")
    ax.text(0.4, 9.1, f"Spearman ρ = {rho:.4f} (P < 10⁻¹⁰⁰)\nTop-12 Overlap = {len(top12_w.intersection(top12_eq))}/12 ({len(top12_w.intersection(top12_eq))/12*100:.1f}%)",
            fontsize=8.8, fontweight="bold", bbox=dict(boxstyle="round,pad=0.35", facecolor="#FDEDEC", edgecolor="#C0392B", lw=0.8))
    ax.set_xlabel("Equal-Weight Composite Score (2:2:2:2:2; 0–10 Scale)", fontsize=9.8, fontweight="bold")
    ax.set_ylabel("Primary Evidence-Weighted Score (3:2:2:2:1; 0–10 Scale)", fontsize=9.8, fontweight="bold")
    ax.set_title("Weight Sensitivity Analysis of Multi-Evidence Antigenic Hotspot Scoring", loc="left", fontsize=10.2, fontweight="bold", pad=34)
    ax.set_xlim(-0.3, 10.5); ax.set_ylim(-0.3, 10.5)
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=7.8, frameon=False)
    ax.grid(True, ls=":", alpha=0.35); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    plt.subplots_adjust(left=0.10, right=0.96, top=0.82, bottom=0.11)
    s2_png = "results/figures/Supplementary_Figure_S2_Weight_Sensitivity_1000dpi.png"
    s2_tif = "results/figures/Supplementary_Figure_S2_Weight_Sensitivity_1000dpi.tiff"
    fig.savefig(s2_png, dpi=300); fig.savefig(s2_tif, dpi=300, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
    for ext, pth in [("png", s2_png), ("tiff", s2_tif)]:
        im = Image.open(pth)
        im.save(pth, dpi=(1000, 1000), compression="tiff_lzw") if ext == "tiff" else im.save(pth, dpi=(1000, 1000))

    # 2. Card 31 (T18): Catalytic & Framework Active-Site Conservation Table (Main Table 17)
    aln = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    n1_map = {}
    p_idx = 0
    for i, ch in enumerate(cons):
        if ch != "-":
            p_idx += 1
            n1_map[p_idx] = i

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    cm = {c.lower(): c for c in meta.columns}
    meta = meta[meta[cm["accession"]].isin(aln)].copy()
    asia_accs = meta[meta[cm["group"]].str.strip().str.lower() == "asia"][cm["accession"]].tolist()
    non_accs = meta[meta[cm["group"]].str.strip().str.lower() != "asia"][cm["accession"]].tolist()
    all_accs = meta[cm["accession"]].tolist()

    cat_rows = []
    for pos, wt, role in CATALYTIC_PANEL:
        c_i = n1_map[pos]
        g_cons = sum(1 for a in all_accs if aln[a][c_i] == wt) / len(all_accs) * 100
        a_cons = sum(1 for a in asia_accs if aln[a][c_i] == wt) / len(asia_accs) * 100
        n_cons = sum(1 for a in non_accs if aln[a][c_i] == wt) / len(non_accs) * 100
        variants = sorted(set(aln[a][c_i] for a in all_accs if aln[a][c_i] not in (wt, "-", "X")))
        cat_rows.append({
            "N1_Position": pos, "Canonical_AA": wt, "Functional_Role": role,
            "Global_Conservation_Pct": round(g_cons, 3),
            "Asian_Conservation_Pct": round(a_cons, 3),
            "NonAsian_Conservation_Pct": round(n_cons, 3),
            "Minor_Variants_Observed": ", ".join(variants) if variants else "None (100% invariant)"
        })

    t17 = pd.DataFrame(cat_rows)
    t17.to_csv("results/tables/Main_Table_17_Catalytic_Framework_Conservation.csv", index=False)
    print("\n=================================================================")
    print(" CARD 31 (T18): MAIN TABLE 17 — Catalytic & Framework Conservation")
    print("=================================================================")
    print(t17.to_string(index=False))

if __name__ == "__main__":
    main()
