#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.patches import Patch

IGNORE_AA = set(["-", "X", "?", ".", "*"])

def read_aligned_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        raise ValueError(f"No sequences found in {fasta_file}")

    seqs = [str(rec.seq).strip().upper() for rec in records]
    lengths = set(len(s) for s in seqs)

    if len(lengths) != 1:
        raise ValueError(f"{fasta_file} is not aligned. All sequences must have equal length.")

    return seqs

def compute_site_variability(seqs, coverage_cutoff=0.80):
    """
    Variation score = 1 - consensus frequency
    Classification:
      conserved     : variation <= 0.02
      variable      : 0.02 < variation < 0.10
      hypervariable : variation >= 0.10
    """
    n_seq = len(seqs)
    aln_len = len(seqs[0])
    ref = seqs[0]

    col_to_refpos = {}
    ref_pos = 0

    for col, aa in enumerate(ref):
        if aa not in IGNORE_AA:
            ref_pos += 1
            col_to_refpos[col] = ref_pos

    rows = []

    for col in range(aln_len):
        if col not in col_to_refpos:
            continue

        residues = []
        for s in seqs:
            aa = s[col]
            if aa not in IGNORE_AA:
                residues.append(aa)

        coverage_frac = len(residues) / n_seq if n_seq > 0 else 0.0
        if coverage_frac < coverage_cutoff:
            continue

        counts = pd.Series(residues).value_counts()
        consensus_aa = counts.idxmax()
        consensus_freq = counts.max() / len(residues)
        variation = 1.0 - consensus_freq

        if variation <= 0.02:
            category = "Conserved"
        elif variation < 0.10:
            category = "Variable"
        else:
            category = "Hypervariable"

        rows.append({
            "Position": col_to_refpos[col],
            "Consensus_AA": consensus_aa,
            "Consensus_frequency": round(consensus_freq, 4),
            "Variation": round(variation, 4),
            "Coverage_pct": round(coverage_frac * 100, 2),
            "Category": category
        })

    return pd.DataFrame(rows)

def analyse_period(label, fasta_file, coverage_cutoff):
    seqs = read_aligned_fasta(fasta_file)
    site_df = compute_site_variability(seqs, coverage_cutoff=coverage_cutoff)

    if site_df.empty:
        raise ValueError(f"No valid analysed positions found for {label}. Check alignment or coverage cutoff.")

    counts = site_df["Category"].value_counts().reindex(
        ["Conserved", "Variable", "Hypervariable"],
        fill_value=0
    )

    total_sites = int(counts.sum())
    percentages = (counts / total_sites * 100).round(1)

    summary = {
        "Period": label,
        "Sequence_count": len(seqs),
        "Analysed_sites": total_sites,
        "Conserved_count": int(counts["Conserved"]),
        "Variable_count": int(counts["Variable"]),
        "Hypervariable_count": int(counts["Hypervariable"]),
        "Conserved_pct": float(percentages["Conserved"]),
        "Variable_pct": float(percentages["Variable"]),
        "Hypervariable_pct": float(percentages["Hypervariable"]),
        "Mean_variation": round(site_df["Variation"].mean(), 3),
        "Max_variation": round(site_df["Variation"].max(), 3)
    }

    safe = label.replace("-", "_")
    site_df.to_csv(f"{safe}_NA469_site_variability_final.csv", index=False)

    return summary, site_df

def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    })

def add_value_labels_stacked(ax, xvals, bottoms, heights, labels, color="white", min_height=4):
    for x, b, h, txt in zip(xvals, bottoms, heights, labels):
        if h >= min_height:
            ax.text(x, b + h/2, txt, ha="center", va="center", color=color, fontsize=10)

def add_bar_labels(ax, bars, fmt="{:.1f}", dy=1.0, fontsize=9):
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            h + dy,
            fmt.format(h),
            ha="center",
            va="bottom",
            fontsize=fontsize
        )

def plot_figure(summary_df, out_prefix):
    set_style()

    periods = summary_df["Period"].tolist()
    x = np.arange(len(periods))

    conserved_pct = summary_df["Conserved_pct"].values
    variable_pct = summary_df["Variable_pct"].values
    hyper_pct = summary_df["Hypervariable_pct"].values

    conserved_count = summary_df["Conserved_count"].values
    variable_count = summary_df["Variable_count"].values
    hyper_count = summary_df["Hypervariable_count"].values

    mean_var = summary_df["Mean_variation"].values
    max_var = summary_df["Max_variation"].values

    xticklabels = [
        f"{p}\n(n={n})"
        for p, n in zip(summary_df["Period"], summary_df["Sequence_count"])
    ]

    colors = {
        "Conserved": "#2F6DB3",
        "Variable": "#BDBDBD",
        "Hypervariable": "#C62828"
    }

    fig = plt.figure(figsize=(16, 12), facecolor="white")
    gs = fig.add_gridspec(3, 2, height_ratios=[1.25, 1.0, 1.0])

    axA = fig.add_subplot(gs[0, :])
    axB = fig.add_subplot(gs[1, 0])
    axC = fig.add_subplot(gs[1, 1])
    axD = fig.add_subplot(gs[2, :])

    # =========================
    # Panel A: Stacked percent
    # =========================
    widthA = 0.8
    axA.bar(x, conserved_pct, width=widthA, color=colors["Conserved"], edgecolor="black", linewidth=0.7, label="Conserved site")
    axA.bar(x, variable_pct, bottom=conserved_pct, width=widthA, color=colors["Variable"], edgecolor="black", linewidth=0.7, label="Variable site")
    axA.bar(x, hyper_pct, bottom=conserved_pct + variable_pct, width=widthA, color=colors["Hypervariable"], edgecolor="black", linewidth=0.7, label="Hypervariable site")

    axA.set_xticks(x)
    axA.set_xticklabels(xticklabels)
    axA.set_ylabel("Classified sites (%)")
    axA.set_ylim(0, 100)
    axA.set_title("Temporal comparison of conserved, variable and hypervariable sites in Asian H1N1 neuraminidase (NA) protein", pad=14)
    axA.grid(axis="y", linestyle="--", alpha=0.25)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.text(-0.055, 1.03, "A", transform=axA.transAxes, fontsize=16, fontweight="bold")

    add_value_labels_stacked(axA, x, np.zeros_like(conserved_pct), conserved_pct, [f"{v:.1f}%" for v in conserved_pct], color="white", min_height=8)
    add_value_labels_stacked(axA, x, conserved_pct, variable_pct, [f"{v:.1f}%" for v in variable_pct], color="black", min_height=5)
    add_value_labels_stacked(axA, x, conserved_pct + variable_pct, hyper_pct, [f"{v:.1f}%" for v in hyper_pct], color="white", min_height=1)

    # legend on top center
    handlesA = [
        Patch(facecolor=colors["Conserved"], edgecolor="black", label="Conserved site"),
        Patch(facecolor=colors["Variable"], edgecolor="black", label="Variable site"),
        Patch(facecolor=colors["Hypervariable"], edgecolor="black", label="Hypervariable site")
    ]
    axA.legend(handles=handlesA, loc="upper center", bbox_to_anchor=(0.5, 1.32), ncol=3, frameon=False)

    # =========================
    # Panel B: Grouped percent
    # =========================
    w = 0.24
    bars1 = axB.bar(x - w, conserved_pct, width=w, color=colors["Conserved"], edgecolor="black", linewidth=0.7, label="Conserved")
    bars2 = axB.bar(x, variable_pct, width=w, color=colors["Variable"], edgecolor="black", linewidth=0.7, label="Variable")
    bars3 = axB.bar(x + w, hyper_pct, width=w, color=colors["Hypervariable"], edgecolor="black", linewidth=0.7, label="Hypervariable")

    axB.set_xticks(x)
    axB.set_xticklabels(xticklabels)
    axB.set_ylabel("Sites (%)")
    axB.set_title("Percentage comparison of all site categories", pad=10)
    axB.grid(axis="y", linestyle="--", alpha=0.25)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.text(-0.13, 1.03, "B", transform=axB.transAxes, fontsize=16, fontweight="bold")
    axB.legend(frameon=False, loc="upper right")

    add_bar_labels(axB, bars1, fmt="{:.1f}", dy=1.2, fontsize=9)
    add_bar_labels(axB, bars2, fmt="{:.1f}", dy=1.2, fontsize=9)
    add_bar_labels(axB, bars3, fmt="{:.1f}", dy=1.2, fontsize=9)

    # =========================
    # Panel C: Grouped counts
    # =========================
    bc1 = axC.bar(x - w, conserved_count, width=w, color=colors["Conserved"], edgecolor="black", linewidth=0.7, label="Conserved")
    bc2 = axC.bar(x, variable_count, width=w, color=colors["Variable"], edgecolor="black", linewidth=0.7, label="Variable")
    bc3 = axC.bar(x + w, hyper_count, width=w, color=colors["Hypervariable"], edgecolor="black", linewidth=0.7, label="Hypervariable")

    axC.set_xticks(x)
    axC.set_xticklabels(xticklabels)
    axC.set_ylabel("Number of sites")
    axC.set_title("Absolute site counts across time periods", pad=10)
    axC.grid(axis="y", linestyle="--", alpha=0.25)
    axC.spines["top"].set_visible(False)
    axC.spines["right"].set_visible(False)
    axC.text(-0.13, 1.03, "C", transform=axC.transAxes, fontsize=16, fontweight="bold")
    axC.legend(frameon=False, loc="upper right")

    add_bar_labels(axC, bc1, fmt="{:.0f}", dy=5, fontsize=9)
    add_bar_labels(axC, bc2, fmt="{:.0f}", dy=5, fontsize=9)
    add_bar_labels(axC, bc3, fmt="{:.0f}", dy=5, fontsize=9)

    # =========================
    # Panel D: Mean + Max variation
    # =========================
    axD.plot(x, mean_var, marker="o", linewidth=2.0, markersize=7, label="Mean variation score")
    axD.plot(x, max_var, marker="s", linewidth=2.0, markersize=7, label="Maximum variation score")

    axD.set_xticks(x)
    axD.set_xticklabels(xticklabels)
    axD.set_ylabel("Variation score")
    axD.set_xlabel("Sequence collection period")
    axD.set_title("Temporal variation score trends", pad=10)
    axD.grid(axis="y", linestyle="--", alpha=0.25)
    axD.spines["top"].set_visible(False)
    axD.spines["right"].set_visible(False)
    axD.text(-0.055, 1.03, "D", transform=axD.transAxes, fontsize=16, fontweight="bold")
    axD.legend(frameon=False, loc="upper left", ncol=2)

    for xi, yi in zip(x, mean_var):
        axD.text(xi, yi + 0.01, f"{yi:.3f}", ha="center", va="bottom", fontsize=9)
    for xi, yi in zip(x, max_var):
        axD.text(xi, yi + 0.01, f"{yi:.3f}", ha="center", va="bottom", fontsize=9)

    fig.subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.10, hspace=0.38, wspace=0.20)

    fig.text(
        0.07,
        0.03,
        "Site classification: conserved = variation ≤ 0.02; variable = 0.02 < variation < 0.10; hypervariable = variation ≥ 0.10. "
        "Positions with <80% non-gap coverage were excluded. The 2021–2025 group contained fewer sequences and should be interpreted cautiously.",
        fontsize=9,
        style="italic"
    )

    plt.savefig(f"{out_prefix}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight", facecolor="white")
    plt.savefig(f"{out_prefix}.svg", bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"[DONE] Saved figure: {out_prefix}.png")
    print(f"[DONE] Saved figure: {out_prefix}.pdf")
    print(f"[DONE] Saved figure: {out_prefix}.svg")

def main():
    parser = argparse.ArgumentParser(description="Generate a final journal-level temporal comparison figure for NA variability.")
    parser.add_argument("--f1", required=True, help="Aligned FASTA for 2015-2017")
    parser.add_argument("--f2", required=True, help="Aligned FASTA for 2018-2020")
    parser.add_argument("--f3", required=True, help="Aligned FASTA for 2021-2025")
    parser.add_argument("--l1", default="2015-2017")
    parser.add_argument("--l2", default="2018-2020")
    parser.add_argument("--l3", default="2021-2025")
    parser.add_argument("-o", "--out", default="H1N1_NA_469_yearwise_final_journal")
    parser.add_argument("--coverage_cutoff", type=float, default=0.80)

    args = parser.parse_args()

    periods = [(args.l1, args.f1), (args.l2, args.f2), (args.l3, args.f3)]

    summaries = []
    all_site_tables = {}

    for label, fasta_file in periods:
        print(f"[INFO] Analysing {label}: {fasta_file}")
        summary, site_df = analyse_period(label, fasta_file, coverage_cutoff=args.coverage_cutoff)
        summaries.append(summary)
        all_site_tables[label] = site_df

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(f"{args.out}_summary.csv", index=False)

    print("\n[SUMMARY]")
    print(summary_df.to_string(index=False))

    plot_figure(summary_df, args.out)

if __name__ == "__main__":
    main()
