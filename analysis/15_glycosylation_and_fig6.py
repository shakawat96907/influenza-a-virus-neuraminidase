#!/usr/bin/env python3
import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

def main():
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, yr_col = col_map["accession"], col_map["group"], col_map["year"]
    meta = meta[meta[acc_col].isin(aln_recs)].copy()
    meta["Is_Asia"] = meta[grp_col].str.strip().str.lower() == "asia"
    meta["Year_Int"] = pd.to_numeric(meta[yr_col], errors="coerce").astype(int)

    def get_p3(y):
        if y <= 2017: return "2015–2017"
        elif y <= 2020: return "2018–2020"
        return "2021–2025"

    meta["Period_3"] = meta["Year_Int"].apply(get_p3)

    # Map canonical N1 position (1..469) -> alignment column index
    n1_to_aln = {}
    n1 = 0
    for idx, ch in enumerate(cons_gapped):
        if ch != "-":
            n1 += 1
            n1_to_aln[n1] = idx

    # Scan N-X-S/T (X != P) across all canonical N1 positions 1..467
    sequon_hits = {pos: {} for pos in range(1, 468)}
    motif_strings = {pos: Counter() for pos in range(1, 468)}

    for acc in meta[acc_col]:
        seq = aln_recs[acc]
        for pos in range(1, 468):
            a1 = seq[n1_to_aln[pos]]
            a2 = seq[n1_to_aln[pos + 1]]
            a3 = seq[n1_to_aln[pos + 2]]
            has_seq = (a1 == "N") and (a2 not in ("P", "-", "X")) and (a3 in ("S", "T"))
            sequon_hits[pos][acc] = has_seq
            motif_strings[pos][f"{a1}{a2}{a3}"] += 1

    # Filter positions with >= 1.0% overall presence in at least one period or overall
    active_positions = []
    for pos in range(1, 468):
        tot_cnt = sum(sequon_hits[pos].values())
        if (tot_cnt / len(meta)) * 100 >= 1.0:
            active_positions.append(pos)
        else:
            # Also check if any single period or year has >= 5.0%
            for p3 in ["2015–2017", "2018–2020", "2021–2025"]:
                sub = meta[meta["Period_3"] == p3][acc_col]
                if len(sub) > 0 and (sum(sequon_hits[pos][a] for a in sub) / len(sub)) * 100 >= 5.0:
                    active_positions.append(pos)
                    break

    t10_rows = []
    yr_rows = []
    years = sorted(meta["Year_Int"].unique())

    for pos in sorted(active_positions):
        dom_motif = motif_strings[pos].most_common(1)
        alt_motifs = "; ".join(f"{m} ({c})" for m, c in motif_strings[pos].most_common(4))

        p_freqs = {}
        for p3 in ["2015–2017", "2018–2020", "2021–2025"]:
            sub = meta[meta["Period_3"] == p3][acc_col]
            p_freqs[p3] = (sum(sequon_hits[pos][a] for a in sub) / len(sub)) * 100

        # Post-2023 powered check (2023-2025)
        sub_2325 = meta[meta["Year_Int"] >= 2023][acc_col]
        f_2325 = (sum(sequon_hits[pos][a] for a in sub_2325) / len(sub_2325)) * 100

        a_sub = meta[meta["Is_Asia"]][acc_col]
        na_sub = meta[~meta["Is_Asia"]][acc_col]
        f_asia = (sum(sequon_hits[pos][a] for a in a_sub) / len(a_sub)) * 100
        f_non = (sum(sequon_hits[pos][a] for a in na_sub) / len(na_sub)) * 100
        tot_n = sum(sequon_hits[pos].values())

        f1, f2, f3 = p_freqs["2015–2017"], p_freqs["2018–2020"], p_freqs["2021–2025"]

        # Exact reproducible classification thresholds
        if f1 >= 80.0 and f2 >= 80.0 and f3 >= 80.0:
            cat = "Persistent (≥80% across all periods)"
        elif f1 < 10.0 and (f3 >= 50.0 or f_2325 >= 50.0):
            cat = "Gained (<10% in 2015–17 to ≥50% post-2020)"
        elif f1 >= 50.0 and (f3 < 10.0 or f_2325 < 10.0):
            cat = "Lost (≥50% in 2015–17 to <10% post-2020)"
        elif f2 >= 10.0 and f1 < 10.0 and f_2325 < 10.0:
            cat = "Transient (≥10% mid-decade, <10% ends)"
        else:
            cat = "Low-frequency / Variable (<10%)"

        t10_rows.append({
            "N1_Sequon_Pos": pos,
            "Dominant_Motif": dom_motif,
            "Total_Count": tot_n,
            "Global_Pct": round(tot_n / len(meta) * 100, 2),
            "Asian_Pct": round(f_asia, 2),
            "NonAsian_Pct": round(f_non, 2),
            "Pct_2015_2017": round(f1, 2),
            "Pct_2018_2020": round(f2, 2),
            "Pct_2021_2025": round(f3, 2),
            "Pct_2023_2025": round(f_2325, 2),
            "Temporal_Classification": cat,
            "Observed_Triplets": alt_motifs
        })

        for yr in years:
            for coh, ia in [("Asia", True), ("Non-Asia", False)]:
                ysub = meta[(meta["Year_Int"] == yr) & (meta["Is_Asia"] == ia)][acc_col]
                ycnt = sum(sequon_hits[pos][a] for a in ysub)
                yr_rows.append({
                    "N1_Sequon_Pos": pos, "Year": yr, "Cohort": coh,
                    "Total_N": len(ysub), "Sequon_Count": ycnt,
                    "Sequon_Pct": round(ycnt / len(ysub) * 100, 2) if len(ysub) > 0 else 0.0
                })

    t10_df = pd.DataFrame(t10_rows)
    t10_df.to_csv("results/tables/Main_Table_10_Glycosylation_Sequons.csv", index=False)
    yr_df = pd.DataFrame(yr_rows)
    yr_df.to_csv("results/tables/Supplementary Table 12 - (annual_glycosylation_sequons_2015_2025).csv", index=False)

    print("=================================================================")
    print(" MAIN TABLE 10: Predicted N-Linked Glycosylation Sequons (N-X-S/T)")
    print("=================================================================")
    print(t10_df.drop(columns=["Observed_Triplets"]).to_string(index=False))

    print("\n=== Detailed Residue & Triplet Breakdown for Key Dynamic Sites (42, 44, 50, 386) ===")
    for p_check in (42, 44, 50, 70, 386):
        cnt44 = sum(sequon_hits[p_check].values())
        trips = '; '.join(f'{m} ({c})' for m, c in motif_strings[p_check].most_common(5))
        print(f'Pos {p_check}: Sequon_Count={cnt44}/1449 ({cnt44/1449*100:.2f}%) | Triplets: {trips}')

    # Render Clean, Overlap-Free Figure 6 (Legends Outside Axes)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6.8), gridspec_kw={"width_ratios": [1.05, 1.15], "wspace": 0.25})

    # Panel A: Period-by-Period Prevalence Across Detected Sequons
    x = np.arange(len(t10_df))
    w = 0.26
    b1 = ax1.bar(x - w, t10_df["Pct_2015_2017"], width=w, color="#2E86C1", label="2015–2017 (n = 497)")
    b2 = ax1.bar(x, t10_df["Pct_2018_2020"], width=w, color="#F39C12", label="2018–2020 (n = 490)")
    b3 = ax1.bar(x + w, t10_df["Pct_2021_2025"], width=w, color="#C0392B", label="2021–2025 (n = 462)")

    ax1.axhline(80, color="#555555", linestyle="--", lw=0.8)
    ax1.axhline(50, color="#7F8C8D", linestyle=":", lw=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"N{p}" for p in t10_df["N1_Sequon_Pos"]], rotation=35, ha="right", fontsize=9, fontweight="bold")
    ax1.set_ylabel("Sequon Prevalence (%)", fontsize=10, fontweight="bold")
    ax1.set_xlabel("Canonical N1 Sequon Position (N-X-S/T, X ≠ P)", fontsize=10, fontweight="bold")
    ax1.set_title("A   Predicted N-Linked Glycosylation Sequon Prevalence Across Periods", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax1.set_ylim(0, 112)
    ax1.legend(handles=[b1, b2, b3], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=8.2, frameon=False)
    ax1.grid(axis="y", linestyle=":", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Annual Trajectory (2015-2025) of Dynamic Sequons (42, 44, 50, 386)
    dyn_styles = [
        (44, "#C0392B", "o", "-", "N44 (Gained)"),
        (386, "#1F618D", "s", "-", "N386 (Lost)"),
        (42, "#27AE60", "^", "--", "N42 (Transient)"),
        (50, "#8E44AD", "d", "--", "N50 (Transient)")
    ]
    handles_b = []
    for pos_d, col_d, mk_d, ls_d, lbl_d in dyn_styles:
        if pos_d in yr_df["N1_Sequon_Pos"].values:
            sub_y = yr_df[yr_df["N1_Sequon_Pos"] == pos_d].groupby("Year").apply(
                lambda g: (g["Sequon_Count"].sum() / g["Total_N"].sum()) * 100
            )
            ln, = ax2.plot(sub_y.index, sub_y.values, color=col_d, marker=mk_d, linestyle=ls_d, lw=2.0, markersize=5.5, label=lbl_d)
            handles_b.append(ln)

    ax2.axvspan(2020.6, 2022.4, color="#F9E79F", alpha=0.4, label="2021–2022 Bottleneck")
    ax2.set_xticks(years)
    ax2.set_xlabel("Collection Year (2015–2025)", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Global Sequon Prevalence (%)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Annual Dynamics of Gained, Lost, and Transient Sequons (2015–2025)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax2.set_ylim(-4, 108)
    ax2.legend(handles=handles_b, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=4, fontsize=8.2, frameon=False)
    ax2.grid(True, linestyle=":", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.97, top=0.87, bottom=0.13)
    fig.savefig("results/figures/Figure6_Glycosylation_Sequons_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure6_Glycosylation_Sequons_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved New Figure 6 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
