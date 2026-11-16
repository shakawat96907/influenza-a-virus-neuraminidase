#!/usr/bin/env python3
import os
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

def compute_variant_metrics(seq_list, cons_seq):
    n = len(seq_list)
    if n == 0:
        return 0, 0, 0, 0.0
    total_vars, rec2_vars, pct1_vars = 0, 0, 0
    total_diffs = 0
    for idx, c_aa in enumerate(cons_seq):
        if c_aa == "-":
            continue
        vars_at_site = [s[idx] for s in seq_list if s[idx] not in ("-", "X", c_aa)]
        cnt = len(vars_at_site)
        total_diffs += cnt
        if cnt >= 1:
            total_vars += 1
        if cnt >= 2:
            rec2_vars += 1
        if (cnt / n) >= 0.01:
            pct1_vars += 1
    mean_diffs_per_seq = total_diffs / n
    return total_vars, rec2_vars, pct1_vars, mean_diffs_per_seq

def rarefy_cohort(seq_list, cons_seq, target_n, n_iter=1000, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    n = len(seq_list)
    if n < target_n:
        t, r2, p1, md = compute_variant_metrics(seq_list, cons_seq)
        return (t, t, t), (r2, r2, r2), (p1, p1, p1), (md, md, md)
    t_arr, r2_arr, p1_arr, md_arr = [], [], [], []
    for _ in range(n_iter):
        idx = rng.choice(n, size=target_n, replace=False)
        sub = [seq_list[i] for i in idx]
        t, r2, p1, md = compute_variant_metrics(sub, cons_seq)
        t_arr.append(t)
        r2_arr.append(r2)
        p1_arr.append(p1)
        md_arr.append(md)
    def summ(a):
        return float(np.mean(a)), float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))
    return summ(t_arr), summ(r2_arr), summ(p1_arr), summ(md_arr)

def main():
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)
    rng = np.random.default_rng(42)

    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons_seq = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, yr_col = col_map["accession"], col_map["group"], col_map["year"]
    meta = meta[meta[acc_col].isin(aln_recs)].copy()
    meta["Is_Asia"] = meta[grp_col].str.strip().str.lower() == "asia"
    meta["Year_Int"] = pd.to_numeric(meta[yr_col], errors="coerce").astype(int)

    def assign_3p(y):
        if y <= 2017: return "2015–2017"
        elif y <= 2020: return "2018–2020"
        return "2021–2025"

    def assign_4p(y):
        if y <= 2017: return "2015–2017 (Pre-clade shift)"
        elif y <= 2020: return "2018–2020 (Pre-pandemic)"
        elif y <= 2022: return "2021–2022 (COVID-19 bottleneck)"
        return "2023–2025 (Post-pandemic surge)"

    meta["Period_3"] = meta["Year_Int"].apply(assign_3p)
    meta["Period_4"] = meta["Year_Int"].apply(assign_4p)

    # 1. Three-Period Comparison (Rarefied to smallest 3-period size = 231)
    print("[1/3] Running 1,000-iteration rarefaction across 3-period and 4-period schemes...")
    p3_order = ["2015–2017", "2018–2020", "2021–2025"]
    min_3p = min(len(meta[(meta["Period_3"] == p) & (meta["Is_Asia"] == ia)]) for p in p3_order for ia in [True, False])

    p4_order = [
        "2015–2017 (Pre-clade shift)",
        "2018–2020 (Pre-pandemic)",
        "2021–2022 (COVID-19 bottleneck)",
        "2023–2025 (Post-pandemic surge)"
    ]
    min_4p = min(len(meta[(meta["Period_4"] == p) & (meta["Is_Asia"] == ia)]) for p in p4_order for ia in [True, False])

    t9_rows = []
    for scheme, col_p, order_p, r_n in [("3-Period Scheme", "Period_3", p3_order, min_3p), ("4-Period Scheme", "Period_4", p4_order, min_4p)]:
        for p in order_p:
            for reg_lbl, ia in [("Asia", True), ("Non-Asia", False)]:
                seqs = [aln_recs[a] for a in meta[(meta[col_p] == p) & (meta["Is_Asia"] == ia)][acc_col]]
                obs_t, obs_r2, obs_p1, obs_md = compute_variant_metrics(seqs, cons_seq)
                (rt_m, rt_lo, rt_hi), (rr2_m, rr2_lo, rr2_hi), _, (rmd_m, rmd_lo, rmd_hi) = rarefy_cohort(seqs, cons_seq, r_n, 1000, rng)
                note = "Descriptive (small N)" if len(seqs) < 50 else "Powered"
                t9_rows.append({
                    "Scheme": scheme, "Period": p, "Cohort": reg_lbl, "N": len(seqs),
                    "Obs_Sites_ge1": obs_t, "Obs_Sites_ge2": obs_r2, "Obs_Sites_ge1pct": obs_p1,
                    "Rarefied_N": r_n,
                    "Rarefied_Sites_ge1_95CI": f"{rt_m:.1f} ({rt_lo:.1f}–{rt_hi:.1f})",
                    "Rarefied_Sites_ge2_95CI": f"{rr2_m:.1f} ({rr2_lo:.1f}–{rr2_hi:.1f})",
                    "Mean_NonCons_Per_Seq": f"{obs_md:.2f}",
                    "Statistical_Status": note,
                    "rt_m": rt_m, "rt_lo": rt_lo, "rt_hi": rt_hi,
                    "rr2_m": rr2_m, "rr2_lo": rr2_lo, "rr2_hi": rr2_hi
                })

    t9_full = pd.DataFrame(t9_rows)
    t9_clean = t9_full.drop(columns=["rt_m", "rt_lo", "rt_hi", "rr2_m", "rr2_lo", "rr2_hi"])
    t9_clean.to_csv("results/tables/Main_Table_9_Temporal_Rarefaction_Summary.csv", index=False)

    print("\n=================================================================")
    print(" MAIN TABLE 9: Temporal Diversity & 1,000-Iter Rarefaction")
    print("=================================================================")
    print(t9_clean.to_string(index=False))

    # 2. Annual (Year-by-Year 2015-2025) Trajectory & Rarefaction (n=20 for powered years)
    print("\n[2/3] Running annual (2015-2025) trajectory & rarefaction...")
    yr_rows = []
    for yr in sorted(meta["Year_Int"].unique()):
        for reg_lbl, ia in [("Asia", True), ("Non-Asia", False)]:
            seqs = [aln_recs[a] for a in meta[(meta["Year_Int"] == yr) & (meta["Is_Asia"] == ia)][acc_col]]
            obs_t, obs_r2, obs_p1, obs_md = compute_variant_metrics(seqs, cons_seq)
            r_target = min(20, len(seqs))
            (rt_m, rt_lo, rt_hi), (rr2_m, rr2_lo, rr2_hi), _, _ = rarefy_cohort(seqs, cons_seq, r_target, 1000, rng)
            yr_rows.append({
                "Year": yr, "Cohort": reg_lbl, "N": len(seqs),
                "Obs_Sites_ge1": obs_t, "Obs_Sites_ge2": obs_r2, "Obs_Sites_ge1pct": obs_p1,
                "Mean_Diffs_Per_Seq": round(obs_md, 2),
                "Rarefied_Target_N": r_target,
                "Rarefied_Sites_ge1": round(rt_m, 1), "Rarefied_95CI_Low": round(rt_lo, 1), "Rarefied_95CI_High": round(rt_hi, 1),
                "Note": "COVID-19 low-surveillance year" if len(seqs) < 15 else "Standard surveillance"
            })
    yr_df = pd.DataFrame(yr_rows)
    yr_df.to_csv("results/tables/Supplementary Table 11 - (annual_temporal_rarefaction_2015_2025).csv", index=False)

    # 3. Render Clean, Overlap-Free Figure 5 (Legends Outside Axes)
    print("\n[3/3] Rendering New Figure 5 (1000 DPI TIFF & 300 DPI PNG)...")
    fig = plt.figure(figsize=(14.8, 9.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], width_ratios=[1.0, 1.05], hspace=0.45, wspace=0.25)

    # Panel A: 3-Period Rarefied Variant Richness (n = 231, 1000 iterations)
    ax1 = fig.add_subplot(gs[0, 0])
    p3_sub = t9_full[t9_full["Scheme"] == "3-Period Scheme"]
    x = np.arange(3)
    w = 0.34
    a_p3 = p3_sub[p3_sub["Cohort"] == "Asia"]
    na_p3 = p3_sub[p3_sub["Cohort"] == "Non-Asia"]

    b1 = ax1.bar(x - w/2, a_p3["rt_m"], width=w, color="#C0392B", yerr=[a_p3["rt_m"] - a_p3["rt_lo"], a_p3["rt_hi"] - a_p3["rt_m"]], capsize=4, label="Asia (Rarefied n = 231)")
    b2 = ax1.bar(x + w/2, na_p3["rt_m"], width=w, color="#1F618D", yerr=[na_p3["rt_m"] - na_p3["rt_lo"], na_p3["rt_hi"] - na_p3["rt_m"]], capsize=4, label="Non-Asia (Rarefied n = 231)")
    for bars in [b1, b2]:
        for bar in bars:
            h = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2, h + 6.5, f"{h:.1f}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(["2015–2017", "2018–2020", "2021–2025"], fontsize=9.5, fontweight="bold")
    ax1.set_ylabel("Rarefied Polymorphic Sites (95% CI)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A   Equal-Size Rarefied Variant Richness Across Periods (n = 231)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax1.set_ylim(0, 240)
    ax1.legend(handles=[b1, b2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8.5, frameon=False)
    ax1.grid(axis="y", linestyle=":", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: 4-Period Decomposition isolating 2021-2022 COVID-19 Bottleneck (n = 36)
    ax2 = fig.add_subplot(gs[0, 1])
    p4_sub = t9_full[t9_full["Scheme"] == "4-Period Scheme"]
    x4 = np.arange(4)
    a_p4 = p4_sub[p4_sub["Cohort"] == "Asia"]
    na_p4 = p4_sub[p4_sub["Cohort"] == "Non-Asia"]

    b3 = ax2.bar(x4 - w/2, a_p4["rr2_m"], width=w, color="#E74C3C", yerr=[a_p4["rr2_m"] - a_p4["rr2_lo"], a_p4["rr2_hi"] - a_p4["rr2_m"]], capsize=3.5, label="Asia (Recurrent ≥2, n = 36)")
    b4 = ax2.bar(x4 + w/2, na_p4["rr2_m"], width=w, color="#2980B9", yerr=[na_p4["rr2_m"] - na_p4["rr2_lo"], na_p4["rr2_hi"] - na_p4["rr2_m"]], capsize=3.5, label="Non-Asia (Recurrent ≥2, n = 36)")
    for bars in [b3, b4]:
        for bar in bars:
            h = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2, h + 2.5, f"{h:.1f}", ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax2.set_xticks(x4)
    ax2.set_xticklabels(["2015–17", "2018–20", "2021–22*\n(Bottleneck)", "2023–25"], fontsize=9, fontweight="bold")
    ax2.set_ylabel("Rarefied Recurrent Sites ≥2 (95% CI)", fontsize=9.5, fontweight="bold")
    ax2.set_title("B   Four-Period Decomposition Isolating 2021–2022 Bottleneck (n = 36)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax2.set_ylim(0, 75)
    ax2.legend(handles=[b3, b4], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8.2, frameon=False)
    ax2.grid(axis="y", linestyle=":", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Panel C (Bottom full width): Annual Accumulation of Non-Consensus Substitutions per Sequence (2015-2025)
    ax3 = fig.add_subplot(gs[1, :])
    ya = yr_df[yr_df["Cohort"] == "Asia"]
    yna = yr_df[yr_df["Cohort"] == "Non-Asia"]
    l1, = ax3.plot(ya["Year"], ya["Mean_Diffs_Per_Seq"], marker="o", lw=2.0, color="#C0392B", label="Asian Cohort Mean Substitutions/Seq")
    l2, = ax3.plot(yna["Year"], yna["Mean_Diffs_Per_Seq"], marker="s", lw=2.0, linestyle="--", color="#1F618D", label="Non-Asian Control Mean Substitutions/Seq")
    ax3.axvspan(2020.6, 2022.4, color="#F9E79F", alpha=0.45, label="COVID-19 Low-Surveillance Bottleneck (2021–2022)")

    ax3.set_xticks(ya["Year"])
    ax3.set_xlabel("Collection Year (2015–2025)", fontsize=10, fontweight="bold")
    ax3.set_ylabel("Mean Non-Consensus AA Differences / Seq", fontsize=9.5, fontweight="bold")
    ax3.set_title("C   Annual Evolutionary Divergence Trajectory Relative to Consensus (2015–2025, N = 1,449)", loc="left", fontsize=11, fontweight="bold", pad=26)
    ax3.legend(handles=[l1, l2, plt.Rectangle((0,0),1,1,fc="#F9E79F",alpha=0.5)], labels=["Asian Cohort", "Non-Asian Control", "2021–2022 Bottleneck (Descriptive)"], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=3, fontsize=8.5, frameon=False)
    ax3.grid(True, linestyle=":", alpha=0.4)

    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.96, top=0.91, bottom=0.08)
    fig.savefig("results/figures/Figure5_Temporal_Rarefaction_And_Trajectory_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure5_Temporal_Rarefaction_And_Trajectory_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved New Figure 5 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
