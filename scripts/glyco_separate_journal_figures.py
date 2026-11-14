#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import SeqIO

IGNORE_AA = set(["-", "X", "?", ".", "*"])
GLYCO_THIRD = set(["S", "T"])

def read_aligned_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        raise ValueError(f"No sequences found in {fasta_file}")

    seqs = [str(rec.seq).strip().upper() for rec in records]
    lengths = set(len(s) for s in seqs)

    if len(lengths) != 1:
        raise ValueError(f"{fasta_file} is not aligned. All sequences must have equal length.")

    return seqs

def is_glyco_motif(triplet):
    if len(triplet) != 3:
        return False

    a = triplet[0]
    b = triplet[1]
    c = triplet[2]

    return a == "N" and b != "P" and c in GLYCO_THIRD

def compute_glyco_prevalence(seqs, min_triplet_coverage=0.80):
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

    for col in range(aln_len - 2):
        if col not in col_to_refpos:
            continue

        start_pos = col_to_refpos[col]
        covered_n = 0
        motif_n = 0

        for s in seqs:
            triplet = s[col:col + 3]

            if len(triplet) < 3:
                continue

            if any(ch in IGNORE_AA for ch in triplet):
                continue

            covered_n += 1

            if is_glyco_motif(triplet):
                motif_n += 1

        coverage_frac = covered_n / n_seq

        if coverage_frac < min_triplet_coverage:
            continue

        motif_pct = (motif_n / covered_n) * 100 if covered_n > 0 else 0

        rows.append({
            "Position": start_pos,
            "Covered_n": covered_n,
            "Coverage_pct": round(coverage_frac * 100, 2),
            "Motif_n": motif_n,
            "Motif_pct": round(motif_pct, 2)
        })

    return pd.DataFrame(rows)

def analyse_period(label, fasta_file, min_triplet_coverage):
    seqs = read_aligned_fasta(fasta_file)

    df = compute_glyco_prevalence(
        seqs,
        min_triplet_coverage=min_triplet_coverage
    )

    safe_label = label.replace("-", "_")
    df.to_csv(f"{safe_label}_glyco_prevalence.csv", index=False)

    summary = {
        "Period": label,
        "Sequence_count": len(seqs),
        "Alignment_length": len(seqs[0]),
        "Analysed_positions": len(df),
        "Mean_motif_pct": round(df["Motif_pct"].mean(), 2),
        "Max_motif_pct": round(df["Motif_pct"].max(), 2)
    }

    return summary, df

def classify_trend(values, threshold):
    present = [v >= threshold for v in values]

    if all(present):
        return "Persistent"
    if (not present[0]) and present[-1]:
        return "Gain"
    if present[0] and (not present[-1]):
        return "Loss"

    return "Transient"

def build_combined(period_labels, period_dfs, presence_threshold, min_display_threshold):
    all_positions = sorted(set().union(*[set(df["Position"]) for df in period_dfs.values()]))

    combined = pd.DataFrame({"Position": all_positions})

    for label in period_labels:
        safe_label = label.replace("-", "_")
        temp = period_dfs[label][["Position", "Motif_pct", "Motif_n", "Covered_n"]].copy()
        temp.columns = [
            "Position",
            f"{safe_label}_pct",
            f"{safe_label}_motif_n",
            f"{safe_label}_covered_n"
        ]
        combined = combined.merge(temp, on="Position", how="left")

    combined = combined.fillna(0)

    pct_cols = [f"{label.replace('-', '_')}_pct" for label in period_labels]

    combined["Max_pct"] = combined[pct_cols].max(axis=1)
    combined["Min_pct"] = combined[pct_cols].min(axis=1)
    combined["Dynamic_range"] = combined["Max_pct"] - combined["Min_pct"]

    combined["Trend"] = combined.apply(
        lambda row: classify_trend(
            [row[pct_cols[0]], row[pct_cols[1]], row[pct_cols[2]]],
            presence_threshold
        ),
        axis=1
    )

    display_df = combined[combined["Max_pct"] >= min_display_threshold].copy()

    if display_df.empty:
        display_df = combined.sort_values("Max_pct", ascending=False).head(10).copy()

    trend_order = {
        "Persistent": 0,
        "Gain": 1,
        "Loss": 2,
        "Transient": 3
    }

    display_df["Trend_order"] = display_df["Trend"].map(trend_order)

    display_df = display_df.sort_values(
        ["Trend_order", "Max_pct", "Dynamic_range"],
        ascending=[True, False, False]
    )

    combined.to_csv("glyco_all_positions.csv", index=False)
    display_df.to_csv("glyco_displayed_positions.csv", index=False)

    return combined, display_df

def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "legend.fontsize": 9,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    })

def save_figure(prefix):
    plt.savefig(f"{prefix}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{prefix}.pdf", bbox_inches="tight", facecolor="white")
    plt.savefig(f"{prefix}.svg", bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"[DONE] Saved: {prefix}.png")
    print(f"[DONE] Saved: {prefix}.pdf")
    print(f"[DONE] Saved: {prefix}.svg")

def plot_heatmap(display_df, summary_df, period_labels, out_prefix):
    set_style()

    safe_labels = [p.replace("-", "_") for p in period_labels]
    pct_cols = [f"{s}_pct" for s in safe_labels]

    heatmap = display_df[pct_cols].values
    site_labels = display_df["Position"].astype(int).astype(str).tolist()

    fig_height = max(5.5, 0.45 * len(site_labels) + 2)
    fig, ax = plt.subplots(figsize=(8.5, fig_height), facecolor="white")

    im = ax.imshow(heatmap, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)

    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cbar.set_label("Predicted glycosylation prevalence (%)")

    xlabels = [
        f"{row['Period']}\n(n={row['Sequence_count']})"
        for _, row in summary_df.iterrows()
    ]

    ax.set_xticks(np.arange(len(period_labels)))
    ax.set_xticklabels(xlabels)

    ax.set_yticks(np.arange(len(site_labels)))
    ax.set_yticklabels(site_labels)

    ax.set_xlabel("Sequence collection period")
    ax.set_ylabel("NA protein position (motif start)")
    ax.set_title("Prevalence heatmap of potential N-linked glycosylation motifs")

    for i in range(heatmap.shape[0]):
        for j in range(heatmap.shape[1]):
            val = heatmap[i, j]
            text_color = "white" if val >= 70 else "black"
            ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=8, color=text_color)

    fig.text(
        0.02,
        0.01,
        "Potential N-linked glycosylation motif: N-X-S/T, where X is not proline. Values represent percentage of sequences carrying the motif.",
        fontsize=9,
        style="italic"
    )

    save_figure(f"{out_prefix}_Figure1_heatmap")

def plot_trend_counts(display_df, out_prefix):
    set_style()

    trend_order = ["Persistent", "Gain", "Loss", "Transient"]
    trend_counts = display_df["Trend"].value_counts().reindex(trend_order, fill_value=0)

    colors = {
        "Persistent": "#2F6DB3",
        "Gain": "#2E8B57",
        "Loss": "#C62828",
        "Transient": "#BFBFBF"
    }

    fig, ax = plt.subplots(figsize=(7, 5.5), facecolor="white")

    bars = ax.bar(
        trend_counts.index,
        trend_counts.values,
        color=[colors[t] for t in trend_counts.index],
        edgecolor="black",
        linewidth=0.7
    )

    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + 0.15,
            str(int(h)),
            ha="center",
            va="bottom",
            fontsize=11
        )

    ax.set_ylabel("Number of glycosylation sites")
    ax.set_xlabel("Glycosylation evolution category")
    ax.set_title("Gain/loss classification of potential glycosylation sites")
    ax.set_ylim(0, max(trend_counts.values) + 2)

    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        Patch(facecolor=colors["Persistent"], edgecolor="black", label="Persistent: present in all periods"),
        Patch(facecolor=colors["Gain"], edgecolor="black", label="Gain: absent/low initially, present later"),
        Patch(facecolor=colors["Loss"], edgecolor="black", label="Loss: present initially, absent/low later"),
        Patch(facecolor=colors["Transient"], edgecolor="black", label="Transient: inconsistent pattern")
    ]

    ax.legend(handles=legend_handles, frameon=False, loc="upper right", fontsize=8)

    fig.text(
        0.02,
        0.01,
        "Trend classification used a 5% prevalence threshold for motif presence.",
        fontsize=9,
        style="italic"
    )

    save_figure(f"{out_prefix}_Figure2_trend_counts")

def plot_dynamic_trajectories(display_df, summary_df, period_labels, out_prefix, top_n):
    set_style()

    safe_labels = [p.replace("-", "_") for p in period_labels]
    pct_cols = [f"{s}_pct" for s in safe_labels]

    top_df = display_df.sort_values(
        ["Dynamic_range", "Max_pct"],
        ascending=False
    ).head(top_n)

    x = np.arange(len(period_labels))

    fig, ax = plt.subplots(figsize=(9, 6), facecolor="white")

    for _, row in top_df.iterrows():
        pos = int(row["Position"])
        y = [row[pct_cols[0]], row[pct_cols[1]], row[pct_cols[2]]]

        ax.plot(
            x,
            y,
            marker="o",
            linewidth=2,
            markersize=6,
            label=f"Site {pos}"
        )

    xlabels = [
        f"{row['Period']}\n(n={row['Sequence_count']})"
        for _, row in summary_df.iterrows()
    ]

    ax.set_xticks(x)
    ax.set_xticklabels(xlabels)

    ax.set_ylabel("Predicted glycosylation prevalence (%)")
    ax.set_xlabel("Sequence collection period")
    ax.set_title(f"Top {len(top_df)} dynamic potential glycosylation sites")
    ax.set_ylim(0, 105)

    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(frameon=False, ncol=2, loc="best")

    fig.text(
        0.02,
        0.01,
        "Dynamic sites were selected based on the largest prevalence change across the three periods.",
        fontsize=9,
        style="italic"
    )

    save_figure(f"{out_prefix}_Figure3_dynamic_trajectories")

def main():
    parser = argparse.ArgumentParser(
        description="Generate separate journal-level glycosylation gain/loss figures."
    )

    parser.add_argument("--f1", required=True, help="Aligned FASTA for 2015-2017")
    parser.add_argument("--f2", required=True, help="Aligned FASTA for 2018-2020")
    parser.add_argument("--f3", required=True, help="Aligned FASTA for 2021-2025")

    parser.add_argument("--l1", default="2015-2017")
    parser.add_argument("--l2", default="2018-2020")
    parser.add_argument("--l3", default="2021-2025")

    parser.add_argument("-o", "--out", default="H1N1_NA_469_glyco")
    parser.add_argument("--presence_threshold", type=float, default=5.0)
    parser.add_argument("--min_display_threshold", type=float, default=5.0)
    parser.add_argument("--min_triplet_coverage", type=float, default=0.80)
    parser.add_argument("--top_dynamic", type=int, default=8)

    args = parser.parse_args()

    period_labels = [args.l1, args.l2, args.l3]
    period_files = [args.f1, args.f2, args.f3]

    summaries = []
    period_dfs = {}

    for label, fasta_file in zip(period_labels, period_files):
        print(f"[INFO] Analysing {label}: {fasta_file}")

        summary, df = analyse_period(
            label,
            fasta_file,
            min_triplet_coverage=args.min_triplet_coverage
        )

        summaries.append(summary)
        period_dfs[label] = df

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(f"{args.out}_period_summary.csv", index=False)

    combined_df, display_df = build_combined(
        period_labels,
        period_dfs,
        presence_threshold=args.presence_threshold,
        min_display_threshold=args.min_display_threshold
    )

    display_df.to_csv(f"{args.out}_displayed_glyco_sites.csv", index=False)

    trend_counts = display_df["Trend"].value_counts().reindex(
        ["Persistent", "Gain", "Loss", "Transient"],
        fill_value=0
    )
    trend_counts.to_csv(f"{args.out}_trend_counts.csv")

    print("\n[PERIOD SUMMARY]")
    print(summary_df.to_string(index=False))

    print("\n[DISPLAYED GLYCOSYLATION SITES]")
    print(display_df[["Position", "Trend", "Max_pct", "Dynamic_range"]].to_string(index=False))

    plot_heatmap(display_df, summary_df, period_labels, args.out)
    plot_trend_counts(display_df, args.out)
    plot_dynamic_trajectories(display_df, summary_df, period_labels, args.out, args.top_dynamic)

if __name__ == "__main__":
    main()
