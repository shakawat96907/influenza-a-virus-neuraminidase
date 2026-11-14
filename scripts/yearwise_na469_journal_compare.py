#!/usr/bin/env python3

import argparse
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO

IGNORE_AA = set(["-", "X", "?", ".", "*"])

def read_aligned_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        raise ValueError(f"No sequences found in {fasta_file}")

    names = []
    seqs = []

    for rec in records:
        names.append(rec.id)
        seqs.append(str(rec.seq).strip().upper())

    lengths = set(len(s) for s in seqs)
    if len(lengths) != 1:
        raise ValueError(f"Sequences in {fasta_file} are not equal length. Use aligned FASTA only.")

    return names, seqs

def compute_variability(seqs, conserved_cutoff=0.02, hyper_cutoff=0.10, min_non_gap_frac=0.80):
    n_seq = len(seqs)
    aln_len = len(seqs[0])
    ref = seqs[0]

    rows = []
    ref_pos = 0

    for col in range(aln_len):
        ref_aa = ref[col]

        if ref_aa in IGNORE_AA:
            continue

        ref_pos += 1

        column = [s[col] for s in seqs if s[col] not in IGNORE_AA]
        non_gap_count = len(column)
        non_gap_frac = non_gap_count / n_seq

        if non_gap_frac < min_non_gap_frac:
            continue

        counts = Counter(column)
        consensus_aa, consensus_count = counts.most_common(1)[0]

        conservation = consensus_count / non_gap_count
        variation = 1.0 - conservation

        probs = np.array([c / non_gap_count for c in counts.values()])
        entropy = -np.sum(probs * np.log2(probs))

        alt_counts = {aa: c for aa, c in counts.items() if aa != ref_aa}

        if alt_counts:
            top_alt, top_alt_count = sorted(alt_counts.items(), key=lambda x: x[1], reverse=True)[0]
            top_alt_frequency = top_alt_count / non_gap_count
            mutation_label = f"{ref_aa}{ref_pos}{top_alt}"
        else:
            top_alt = ""
            top_alt_count = 0
            top_alt_frequency = 0
            mutation_label = ""

        if variation <= conserved_cutoff:
            category = "Conserved"
        elif variation >= hyper_cutoff:
            category = "Hypervariable"
        else:
            category = "Variable"

        rows.append({
            "Position": ref_pos,
            "Reference_AA": ref_aa,
            "Consensus_AA": consensus_aa,
            "Conservation": round(conservation, 4),
            "Variation": round(variation, 4),
            "Entropy": round(float(entropy), 4),
            "Non_gap_count": non_gap_count,
            "Top_alt_AA": top_alt,
            "Top_alt_count": top_alt_count,
            "Top_alt_frequency": round(top_alt_frequency, 4),
            "Mutation_label": mutation_label,
            "Category": category
        })

    return pd.DataFrame(rows)

def analyze_period(label, fasta_file, conserved_cutoff, hyper_cutoff, min_non_gap_frac):
    names, seqs = read_aligned_fasta(fasta_file)

    df = compute_variability(
        seqs,
        conserved_cutoff=conserved_cutoff,
        hyper_cutoff=hyper_cutoff,
        min_non_gap_frac=min_non_gap_frac
    )

    total_sites = len(df)
    counts = df["Category"].value_counts().to_dict()

    conserved_n = counts.get("Conserved", 0)
    variable_n = counts.get("Variable", 0)
    hyper_n = counts.get("Hypervariable", 0)

    summary = {
        "Period": label,
        "File": fasta_file,
        "Sequence_count": len(seqs),
        "Alignment_length": len(seqs[0]),
        "Classified_sites": total_sites,
        "Conserved_n": conserved_n,
        "Variable_n": variable_n,
        "Hypervariable_n": hyper_n,
        "Conserved_pct": round((conserved_n / total_sites) * 100, 2),
        "Variable_pct": round((variable_n / total_sites) * 100, 2),
        "Hypervariable_pct": round((hyper_n / total_sites) * 100, 2),
        "Mean_variation": round(df["Variation"].mean(), 4),
        "Median_variation": round(df["Variation"].median(), 4),
        "Max_variation": round(df["Variation"].max(), 4)
    }

    safe_label = label.replace("-", "_")
    df.to_csv(f"{safe_label}_NA469_site_variability.csv", index=False)

    top_sites = df.sort_values(["Variation", "Top_alt_frequency"], ascending=False).head(15)
    top_sites.to_csv(f"{safe_label}_NA469_top_variable_sites.csv", index=False)

    return summary, df

def plot_journal_figure(summary_df, all_site_dfs, out_prefix, title):
    colors = {
        "Conserved": "#2F6DB3",
        "Variable": "#BFBFBF",
        "Hypervariable": "#C62828"
    }

    periods = summary_df["Period"].tolist()
    n_labels = [f"{p}\n(n={n})" for p, n in zip(summary_df["Period"], summary_df["Sequence_count"])]

    conserved_pct = summary_df["Conserved_pct"].values
    variable_pct = summary_df["Variable_pct"].values
    hyper_pct = summary_df["Hypervariable_pct"].values

    conserved_n = summary_df["Conserved_n"].values
    variable_n = summary_df["Variable_n"].values
    hyper_n = summary_df["Hypervariable_n"].values

    mean_var = summary_df["Mean_variation"].values
    max_var = summary_df["Max_variation"].values

    x = np.arange(len(periods))
    width = 0.35

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    })

    fig = plt.figure(figsize=(15, 10), facecolor="white")

    # Panel A: stacked percentage bar
    ax1 = fig.add_axes([0.08, 0.62, 0.88, 0.26])

    ax1.bar(x, conserved_pct, color=colors["Conserved"], edgecolor="black", linewidth=0.6)
    ax1.bar(x, variable_pct, bottom=conserved_pct, color=colors["Variable"], edgecolor="black", linewidth=0.6)
    ax1.bar(x, hyper_pct, bottom=conserved_pct + variable_pct, color=colors["Hypervariable"], edgecolor="black", linewidth=0.6)

    for i in range(len(x)):
        ax1.text(x[i], conserved_pct[i] / 2, f"{conserved_pct[i]:.1f}%", ha="center", va="center", fontsize=10, color="white")
        if variable_pct[i] >= 3:
            ax1.text(x[i], conserved_pct[i] + variable_pct[i] / 2, f"{variable_pct[i]:.1f}%", ha="center", va="center", fontsize=9, color="black")
        if hyper_pct[i] >= 1:
            ax1.text(x[i], conserved_pct[i] + variable_pct[i] + hyper_pct[i] / 2, f"{hyper_pct[i]:.1f}%", ha="center", va="center", fontsize=9, color="white")

    ax1.set_xticks(x)
    ax1.set_xticklabels(n_labels)
    ax1.set_ylim(0, 100)
    ax1.set_ylabel("Classified sites (%)")
    ax1.set_title(title, pad=18)
    ax1.grid(axis="y", linestyle="--", alpha=0.25)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.text(-0.055, 1.04, "A", transform=ax1.transAxes, fontsize=15, fontweight="bold")

    # Panel B: conserved vs hypervariable percentage
    ax2 = fig.add_axes([0.08, 0.35, 0.40, 0.20])

    b1 = ax2.bar(x - width / 2, conserved_pct, width, color=colors["Conserved"], edgecolor="black", linewidth=0.6, label="Conserved")
    b2 = ax2.bar(x + width / 2, hyper_pct, width, color=colors["Hypervariable"], edgecolor="black", linewidth=0.6, label="Hypervariable")

    for bar in list(b1) + list(b2):
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2, h + 0.7, f"{h:.1f}%", ha="center", va="bottom", fontsize=9)

    ax2.set_xticks(x)
    ax2.set_xticklabels(n_labels)
    ax2.set_ylabel("Sites (%)")
    ax2.set_ylim(0, max(conserved_pct.max(), hyper_pct.max()) + 8)
    ax2.grid(axis="y", linestyle="--", alpha=0.25)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(-0.13, 1.04, "B", transform=ax2.transAxes, fontsize=15, fontweight="bold")

    # Panel C: absolute site counts
    ax3 = fig.add_axes([0.56, 0.35, 0.40, 0.20])

    ax3.bar(x - width, conserved_n, width, color=colors["Conserved"], edgecolor="black", linewidth=0.6)
    ax3.bar(x, variable_n, width, color=colors["Variable"], edgecolor="black", linewidth=0.6)
    ax3.bar(x + width, hyper_n, width, color=colors["Hypervariable"], edgecolor="black", linewidth=0.6)

    for i in range(len(x)):
        ax3.text(x[i] - width, conserved_n[i] + 2, str(conserved_n[i]), ha="center", va="bottom", fontsize=9)
        ax3.text(x[i], variable_n[i] + 2, str(variable_n[i]), ha="center", va="bottom", fontsize=9)
        ax3.text(x[i] + width, hyper_n[i] + 2, str(hyper_n[i]), ha="center", va="bottom", fontsize=9)

    ax3.set_xticks(x)
    ax3.set_xticklabels(n_labels)
    ax3.set_ylabel("Number of sites")
    ax3.set_ylim(0, max(conserved_n.max(), variable_n.max(), hyper_n.max()) + 30)
    ax3.grid(axis="y", linestyle="--", alpha=0.25)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.text(-0.13, 1.04, "C", transform=ax3.transAxes, fontsize=15, fontweight="bold")

    # Panel D: mean and maximum variation scores
    ax4 = fig.add_axes([0.08, 0.12, 0.88, 0.16])

    ax4.plot(x, mean_var, marker="o", linewidth=2, label="Mean variation score")
    ax4.plot(x, max_var, marker="s", linewidth=2, label="Maximum variation score")

    for i in range(len(x)):
        ax4.text(x[i], mean_var[i] + 0.01, f"{mean_var[i]:.3f}", ha="center", va="bottom", fontsize=9)
        ax4.text(x[i], max_var[i] + 0.01, f"{max_var[i]:.3f}", ha="center", va="bottom", fontsize=9)

    ax4.set_xticks(x)
    ax4.set_xticklabels(n_labels)
    ax4.set_ylabel("Variation score")
    ax4.set_xlabel("Sequence collection period")
    ax4.set_ylim(0, max(max_var.max() + 0.08, 0.25))
    ax4.grid(axis="y", linestyle="--", alpha=0.25)
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    ax4.legend(frameon=False, loc="upper left", ncol=2)
    ax4.text(-0.055, 1.06, "D", transform=ax4.transAxes, fontsize=15, fontweight="bold")

    legend_handles = [
        Patch(facecolor=colors["Conserved"], edgecolor="black", label="Conserved site"),
        Patch(facecolor=colors["Variable"], edgecolor="black", label="Variable site"),
        Patch(facecolor=colors["Hypervariable"], edgecolor="black", label="Hypervariable site")
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=3,
        frameon=False
    )

    fig.text(
        0.08,
        0.035,
        "Site classification: conserved = variation ≤ 0.02; variable = 0.02 < variation < 0.10; hypervariable = variation ≥ 0.10. "
        "Positions with <80% non-gap coverage were excluded. The 2021–2025 group has a smaller sample size and should be interpreted cautiously.",
        fontsize=9,
        style="italic"
    )

    png = f"{out_prefix}.png"
    pdf = f"{out_prefix}.pdf"
    svg = f"{out_prefix}.svg"

    plt.savefig(png, dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.savefig(svg, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"[DONE] Saved figure: {png}")
    print(f"[DONE] Saved figure: {pdf}")
    print(f"[DONE] Saved figure: {svg}")

def main():
    parser = argparse.ArgumentParser(description="Journal-level year-wise comparison of conserved, variable and hypervariable sites in H1N1 NA protein.")
    parser.add_argument("--f1", required=True, help="Aligned FASTA for 2015-2017")
    parser.add_argument("--f2", required=True, help="Aligned FASTA for 2018-2020")
    parser.add_argument("--f3", required=True, help="Aligned FASTA for 2021-2025")
    parser.add_argument("--l1", default="2015-2017")
    parser.add_argument("--l2", default="2018-2020")
    parser.add_argument("--l3", default="2021-2025")
    parser.add_argument("-o", "--out", default="H1N1_NA_469_yearwise_journal")
    parser.add_argument("--title", default="Temporal comparison of conserved and hypervariable sites in Asian H1N1 neuraminidase (NA) protein")
    parser.add_argument("--conserved_cutoff", type=float, default=0.02)
    parser.add_argument("--hyper_cutoff", type=float, default=0.10)
    parser.add_argument("--min_non_gap_frac", type=float, default=0.80)

    args = parser.parse_args()

    summaries = []
    all_dfs = {}

    for label, fasta in [(args.l1, args.f1), (args.l2, args.f2), (args.l3, args.f3)]:
        print(f"\n[INFO] Analysing {label}: {fasta}")
        summary, df = analyze_period(
            label,
            fasta,
            args.conserved_cutoff,
            args.hyper_cutoff,
            args.min_non_gap_frac
        )
        summaries.append(summary)
        all_dfs[label] = df

    summary_df = pd.DataFrame(summaries)
    summary_csv = f"{args.out}_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    print("\n[SUMMARY TABLE]")
    print(summary_df.to_string(index=False))
    print(f"\n[DONE] Saved summary table: {summary_csv}")

    plot_journal_figure(
        summary_df,
        all_dfs,
        args.out,
        args.title
    )

if __name__ == "__main__":
    main()
