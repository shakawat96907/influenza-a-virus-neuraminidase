#!/usr/bin/env python3
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

def parse_hyphy_table(json_path):
    with open(json_path) as f:
        data = json.load(f)
    headers = [next(iter(h)) for h in data["MLE"]["headers"]]
    content = data["MLE"]["content"]["0"]
    return pd.DataFrame([r[:len(headers)] for r in content], columns=headers)

def main():
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    # Map codon index (1..493) -> canonical N1 position (1..469)
    codon_to_n1 = {}
    n1 = 0
    for idx, ch in enumerate(cons_gapped, start=1):
        if ch != "-":
            n1 += 1
            codon_to_n1[idx] = (n1, ch)

    fel = parse_hyphy_table("results/selection/fel_dedup.json")
    fubar = parse_hyphy_table("results/selection/fubar_dedup.json")
    meme = parse_hyphy_table("results/selection/meme_dedup.json")

    rows = []
    for c_idx in range(1, len(fel) + 1):
        if c_idx not in codon_to_n1:
            continue
        n1_pos, cons_aa = codon_to_n1[c_idx]
        r_fel = fel.iloc[c_idx - 1]
        r_fub = fubar.iloc[c_idx - 1]
        r_mem = meme.iloc[c_idx - 1]

        alpha_fel = float(r_fel["alpha"])
        beta_fel = float(r_fel["beta"])
        p_fel = float(r_fel["p-value"])
        omega_fel = round(beta_fel / alpha_fel, 3) if alpha_fel > 1e-6 else (99.0 if beta_fel > 0 else 0.0)

        alpha_fub = float(r_fub["alpha"])
        beta_fub = float(r_fub["beta"])
        pp_pos = float(r_fub["Prob[alpha<beta]"])
        pp_neg = float(r_fub["Prob[alpha>beta]"])
        omega_fub = round(beta_fub / alpha_fub, 3) if alpha_fub > 1e-6 else 0.0

        p_meme = float(r_mem["p-value"])
        branches_meme = int(r_mem["# branches under selection"])

        # Consensus classification across FEL (p<=0.10), MEME (p<=0.10), FUBAR (PP>=0.90)
        pos_signals = []
        if p_fel <= 0.10 and beta_fel > alpha_fel: pos_signals.append("FEL")
        if p_meme <= 0.10: pos_signals.append("MEME")
        if pp_pos >= 0.90: pos_signals.append("FUBAR")

        neg_signals = []
        if p_fel <= 0.10 and beta_fel < alpha_fel: neg_signals.append("FEL")
        if pp_neg >= 0.90: neg_signals.append("FUBAR")

        if pos_signals:
            sel_class = f"Positive / Diversifying ({'+'.join(pos_signals)})"
        elif neg_signals:
            sel_class = f"Purifying / Negative ({'+'.join(neg_signals)})"
        else:
            sel_class = "Neutral / Constraint-buffered"

        rows.append({
            "N1_Position": n1_pos, "Consensus_AA": cons_aa,
            "FEL_alpha_dS": round(alpha_fel, 4), "FEL_beta_dN": round(beta_fel, 4),
            "FEL_dN_dS": omega_fel, "FEL_P": round(p_fel, 4),
            "FUBAR_dN_dS": omega_fub, "FUBAR_PP_Pos": round(pp_pos, 4), "FUBAR_PP_Neg": round(pp_neg, 4),
            "MEME_P": round(p_meme, 4), "MEME_Branches": branches_meme,
            "Selection_Regime": sel_class
        })

    df = pd.DataFrame(rows)
    df.to_csv("results/tables/Supplementary Table 13 - (codon_selection_FEL_MEME_FUBAR).csv", index=False)

    # Summary counts
    n_pos_any = df["Selection_Regime"].str.startswith("Positive").sum()
    n_neg_any = df["Selection_Regime"].str.startswith("Purifying").sum()
    n_neu_any = df["Selection_Regime"].str.startswith("Neutral").sum()
    print(f"Selection Summary across 469 N1 sites -> Positive/Diversifying: {n_pos_any} | Purifying: {n_neg_any} | Neutral: {n_neu_any}")

    # Cross-check Positive sites + our Candidate Variants
    cand_positions = {74, 51, 452, 416, 389, 222, 66, 19, 200, 81, 77, 339, 469, 453, 50}
    pos_positions = set(df[df["Selection_Regime"].str.startswith("Positive")]["N1_Position"])
    target_positions = sorted(cand_positions | pos_positions)

    t11 = df[df["N1_Position"].isin(target_positions)].copy()
    t11["Is_Candidate_Variant"] = t11["N1_Position"].apply(lambda p: "Yes" if p in cand_positions else "Positive-only")
    t11.to_csv("results/tables/Main_Table_11_Selection_And_Candidates.csv", index=False)

    print("\n=================================================================")
    print(" MAIN TABLE 11: Positive Selection Sites & Candidate Cross-Check")
    print("=================================================================")
    print(t11[["N1_Position", "Consensus_AA", "Is_Candidate_Variant", "FEL_dN_dS", "FEL_P", "FUBAR_dN_dS", "FUBAR_PP_Pos", "FUBAR_PP_Neg", "MEME_P", "MEME_Branches", "Selection_Regime"]].to_string(index=False))

    # Render Clean, Overlap-Free Figure 7 (Legends Outside Axes)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14.2, 8.8), gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.44})

    # Panel A: Site-wise FUBAR dN/dS profile across 1-469 colored by selection regime
    neu_df = df[df["Selection_Regime"].str.startswith("Neutral")]
    pur_df = df[df["Selection_Regime"].str.startswith("Purifying")]
    pos_df = df[df["Selection_Regime"].str.startswith("Positive")]

    b_neu = ax1.scatter(neu_df["N1_Position"], neu_df["FUBAR_dN_dS"], s=14, color="#95A5A6", alpha=0.65, label=f"Neutral (n = {len(neu_df)})")
    b_pur = ax1.scatter(pur_df["N1_Position"], pur_df["FUBAR_dN_dS"], s=18, color="#2471A3", alpha=0.80, label=f"Purifying Selection (n = {len(pur_df)})")
    b_pos = ax1.scatter(pos_df["N1_Position"], pos_df["FUBAR_dN_dS"], s=52, color="#C0392B", edgecolors="#111111", lw=0.7, zorder=4, label=f"Positive / Diversifying Selection (n = {len(pos_df)})")
    h_one = ax1.axhline(1.0, color="#333333", linestyle="--", lw=1.0, label="Neutrality Threshold (ω = 1.0)")

    for i_idx, (_, r) in enumerate(pos_df.iterrows()):
        dy = 0.25 if i_idx % 2 == 0 else -0.32
        ax1.annotate(f"{r['Consensus_AA']}{int(r['N1_Position'])}", xy=(r["N1_Position"], r["FUBAR_dN_dS"]), xytext=(r["N1_Position"], max(0.15, r["FUBAR_dN_dS"] + dy)), fontsize=8, fontweight="bold", ha="center", color="#78281F", arrowprops=dict(arrowstyle="-", lw=0.6, color="#555555"))

    ax1.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Site-Wise ω (dN/dS, FUBAR)", fontsize=10, fontweight="bold")
    ax1.set_title("A   Codon-Level Selection Profile Across Human A(H1N1)pdm09 Neuraminidase (FEL, MEME, FUBAR)", loc="left", fontsize=11, fontweight="bold", pad=28)
    ax1.set_xlim(0, 474)
    ax1.set_ylim(-0.1, max(3.6, df["FUBAR_dN_dS"].max() * 1.15))
    ax1.legend(handles=[b_pos, b_pur, b_neu, h_one], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=4, fontsize=8.2, frameon=False)
    ax1.grid(True, linestyle=":", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Episodic Diversifying Selection (-log10 MEME p-value) & Candidate Cross-Check
    neg_log_p = -np.log10(np.clip(df["MEME_P"], 1e-6, 1.0))
    cand_mask = df["N1_Position"].isin(cand_positions)
    s_other = ax2.scatter(df.loc[~cand_mask, "N1_Position"], neg_log_p[~cand_mask], s=14, color="#BDC3C7", alpha=0.7, label="Other N1 Sites")
    s_cand = ax2.scatter(df.loc[cand_mask, "N1_Position"], neg_log_p[cand_mask], s=46, color="#E67E22", edgecolors="#111111", lw=0.7, zorder=4, label="Key Polymorphic / Regional Candidates (n = 15)")
    h_p10 = ax2.axhline(-np.log10(0.10), color="#C0392B", linestyle="--", lw=1.0, label="MEME Significance Threshold (p = 0.10)")

    for _, r in df[cand_mask & (df["MEME_P"] <= 0.10)].iterrows():
        y_val = -np.log10(max(r["MEME_P"], 1e-6))
        ax2.annotate(f"{r['Consensus_AA']}{int(r['N1_Position'])} (p={r['MEME_P']:.3f})", xy=(r["N1_Position"], y_val), xytext=(r["N1_Position"] - 15, y_val + 0.35), fontsize=8.2, fontweight="bold", color="#922B21", arrowprops=dict(arrowstyle="-", lw=0.65, color="#555555"))

    ax2.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax2.set_ylabel("−log₁₀(MEME P-value)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Episodic Diversifying Selection (MEME) Across Candidate Polymorphic & Enriched Sites", loc="left", fontsize=11, fontweight="bold", pad=28)
    ax2.set_xlim(0, 474)
    ax2.set_ylim(-0.15, max(4.3, neg_log_p.max() * 1.15))
    ax2.legend(handles=[s_cand, s_other, h_p10], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=3, fontsize=8.2, frameon=False)
    ax2.grid(True, linestyle=":", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.97, top=0.91, bottom=0.08)
    fig.savefig("results/figures/Figure7_Codon_Selection_Profile_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure7_Codon_Selection_Profile_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved New Figure 7 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
