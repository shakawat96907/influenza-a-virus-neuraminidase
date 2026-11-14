#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
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
    a, b, c = triplet[0], triplet[1], triplet[2]
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
            triplet = s[col:col+3]

            if len(triplet) < 3:
                continue
            if any(ch in IGNORE_AA for ch in triplet):
                continue

            covered_n += 1
            if is_glyco_motif(triplet):
                motif_n += 1

        coverage_frac = covered_n / n_seq if n_seq > 0 else 0.0
        if coverage_frac < min_triplet_coverage:
            continue

        motif_pct = (motif_n / covered_n) * 100 if covered_n > 0 else 0.0

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
    df = compute_glyco_prevalence(seqs, min_triplet_coverage=min_triplet_coverage)

    safe_label = label.replace("-", "_")
    df.to_csv(f"{safe_label}_glyco_prevalence_combined.csv", index=False)

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
    elif (not present[0]) and present[-1]:
        return "Gain"
    elif present[0] and (not present[-1]):
        return "Loss"
    else:
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
        display_df = combined.sort_values(["Max_pct", "Dynamic_range"], ascending=False).head(10).copy()

    trend_order = {"Persistent": 0, "Gain": 1, "Loss": 2, "Transient": 3}
    display_df["Trend_order"] = display_df["Trend"].map(trend_order)

    display_df = display_df.sort_values(
        ["Trend_order", "Max_pct", "Dynamic_range"],
        ascending=[True, False, False]
    )

    combined.to_csv("glyco_combined_all_positions.csv", index=False)
    display_df.to_csv("glyco_combined_displayed_positions.csv", index=False)

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

def plot_combined_figure(display_df, summary_df, period_labels, out_prefix, presence_threshold, top_dynamic):
    set_style()

    safe_labels = [p.replace("-", "_") for p in period_labels]
    pct_cols = [f"{s}_pct" for s in safe_labels]

    heatmap = display_df[pct_cols].values
    site_labels = display_df["Position"].astype(int).astype(str).tolist()

    top_df = display_df.sort_values(
        ["Dynamic_range", "Max_pct"],
        ascending=False
    ).head(top_dynamic)

    trend_order = ["Persistent", "Gain", "Loss", "Transient"]
    trend_counts = display_df["Trend"].value_counts().reindex(trend_order, fill_value=0)

    # make figure large enough
    fig_height = max(11.5, 0.35 * len(site_labels) + 7.5)
    fig = plt.figure(figsize=(16, fig_height), constrained_layout=True, facecolor="white")
    gs = GridSpec(
        3, 3, figure=fig,
        width_ratios=[1.6, 1.6, 1.05],
        height_ratios=[1.05, 1.05, 1.25]
    )

    # ---------- Panel A: Heatmap ----------
    axA = fig.add_subplot(gs[0:2, 0:2])

    im = axA.imshow(heatmap, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)

    xlabels = [
        f"{row['Period']}\n(n={row['Sequence_count']})"
        for _, row in summary_df.iterrows()
    ]

    axA.set_xticks(np.arange(len(period_labels)))
    axA.set_xticklabels(xlabels)
    axA.set_yticks(np.arange(len(site_labels)))
    axA.set_yticklabels(site_labels)

    axA.set_xlabel("Sequence collection period")
    axA.set_ylabel("NA protein position (motif start)")
    axA.set_title("Prevalence heatmap of potential N-linked glycosylation motifs", pad=12)
    axA.text(-0.08, 1.02, "A", transform=axA.transAxes, fontsize=15, fontweight="bold")

    if len(site_labels) <= 18:
        for i in range(heatmap.shape[0]):
            for j in range(heatmap.shape[1]):
                val = heatmap[i, j]
                txt_color = "white" if val >= 70 else "black"
                axA.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=8.5, color=txt_color)

    cbar = fig.colorbar(im, ax=axA, fraction=0.025, pad=0.02)
    cbar.set_label("Predicted glycosylation prevalence (%)")

    # ---------- Panel B: Trend counts ----------
    axB = fig.add_subplot(gs[0:2, 2])

    trend_colors = {
        "Persistent": "#2F6DB3",
        "Gain": "#2E8B57",
        "Loss": "#C62828",
        "Transient": "#BFBFBF"
    }

    bars = axB.bar(
        trend_counts.index,
        trend_counts.values,
        color=[trend_colors[t] for t in trend_counts.index],
        edgecolor="black",
        linewidth=0.7
    )

    for bar in bars:
        h = bar.get_height()
        axB.text(
            bar.get_x() + bar.get_width() / 2,
            h + 0.12,
            str(int(h)),
            ha="center",
            va="bottom",
            fontsize=10
        )

    axB.set_ylabel("Number of glycosylation sites")
    axB.set_xlabel("Evolution category")
    axB.set_title("Trend category counts", pad=12)
    axB.grid(axis="y", linestyle="--", alpha=0.25)
    axB.spines["top"].set_visible(False)
    axB.spines["right"].set_visible(False)
    axB.set_ylim(0, max(trend_counts.values) + 2)
    axB.tick_params(axis='x', rotation=20)
    axB.text(-0.25, 1.02, "B", transform=axB.transAxes, fontsize=15, fontweight="bold")

    legend_handles = [
        Patch(facecolor=trend_colors["Persistent"], edgecolor="black", label="Persistent"),
        Patch(facecolor=trend_colors["Gain"], edgecolor="black", label="Gain"),
        Patch(facecolor=trend_colors["Loss"], edgecolor="black", label="Loss"),
        Patch(facecolor=trend_colors["Transient"], edgecolor="black", label="Transient")
    ]
    axB.legend(handles=legend_handles, frameon=False, loc="upper right")

    # ---------- Panel C: Dynamic trajectories ----------
    axC = fig.add_subplot(gs[2, 0:3])

    x = np.arange(len(period_labels))

    for _, row in top_df.iterrows():
        pos = int(row["Position"])
        y = [row[pct_cols[0]], row[pct_cols[1]], row[pct_cols[2]]]
        axC.plot(
            x, y,
            marker="o",
            linewidth=2.2,
            markersize=6,
            label=f"Site {pos}"
        )

    axC.set_xticks(x)
    axC.set_xticklabels(xlabels)
    axC.set_ylabel("Predicted glycosylation prevalence (%)")
    axC.set_xlabel("Sequence collection period")
    axC.set_title(f"Top {len(top_df)} dynamic potential glycosylation sites", pad=12)
    axC.set_ylim(0, 105)
    axC.grid(axis="y", linestyle="--", alpha=0.25)
    axC.spines["top"].set_visible(False)
    axC.spines["right"].set_visible(False)
    axC.text(-0.02, 1.02, "C", transform=axC.transAxes, fontsize=15, fontweight="bold")

    # put legend outside so lines/text do not overlap
    axC.legend(
        frameon=False,
        ncol=2,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        title="Dynamic sites"
    )

    fig.suptitle(
        "Evolution map of potential N-linked glycosylation motifs in Asian H1N1 neuraminidase (NA)",
        fontsize=16,
        y=1.02
    )

    fig.text(
        0.01,
        0.005,
        f"Potential N-linked glycosylation motifs were defined as N-X-S/T (X ≠ P). "
        f"Site prevalence represents the percentage of sequences carrying the motif at a given NA position. "
        f"Trend classification used a {presence_threshold:.1f}% prevalence threshold. "
        f"Only sites with maximum prevalence ≥ 5% in at least one period are displayed. "
        f"The 2021–2025 group contains fewer sequences and should be interpreted cautiously.",
        fontsize=9,
        style="italic"
    )

    plt.savefig(f"{out_prefix}.png", dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight", facecolor="white")
    plt.savefig(f"{out_prefix}.svg", bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"[DONE] Saved: {out_prefix}.png")
    print(f"[DONE] Saved: {out_prefix}.pdf")
    print(f"[DONE] Saved: {out_prefix}.svg")

def main():
    parser = argparse.ArgumentParser(description="Generate a large combined journal-level glycosylation evolution figure.")
    parser.add_argument("--f1", required=True, help="Aligned FASTA for 2015-2017")
    parser.add_argument("--f2", required=True, help="Aligned FASTA for 2018-2020")
    parser.add_argument("--f3", required=True, help="Aligned FASTA for 2021-2025")
    parser.add_argument("--l1", default="2015-2017")
    parser.add_argument("--l2", default="2018-2020")
    parser.add_argument("--l3", default="2021-2025")
    parser.add_argument("-o", "--out", default="H1N1_NA_469_glyco_combined_journal")
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
        summary, df = analyse_period(label, fasta_file, args.min_triplet_coverage)
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

    display_df.to_csv(f"{args.out}_displayed_sites.csv", index=False)

    print("\n[PERIOD SUMMARY]")
    print(summary_df.to_string(index=False))

    print("\n[DISPLAYED SITES]")
    print(display_df[["Position", "Trend", "Max_pct", "Dynamic_range"]].to_string(index=False))

    plot_combined_figure(
        display_df=display_df,
        summary_df=summary_df,
        period_labels=period_labels,
        out_prefix=args.out,
        presence_threshold=args.presence_threshold,
        top_dynamic=args.top_dynamic
    )

if __name__ == "__main__":
    main()
