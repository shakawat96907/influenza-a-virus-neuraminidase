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

    names = []
    seqs = []

    for rec in records:
        names.append(rec.id)
        seqs.append(str(rec.seq).strip().upper())

    lengths = set(len(s) for s in seqs)
    if len(lengths) != 1:
        raise ValueError(f"Sequences in {fasta_file} are not aligned / equal length.")

    return names, seqs

def is_n_glyco_motif(triplet):
    """
    N-linked glycosylation sequon rule:
    N-X-S/T where X != P
    """
    if len(triplet) != 3:
        return False
    a, b, c = triplet[0], triplet[1], triplet[2]
    if a == "N" and b != "P" and c in GLYCO_THIRD:
        return True
    return False

def compute_glyco_prevalence(seqs, min_triplet_coverage=0.80):
    n_seq = len(seqs)
    aln_len = len(seqs[0])
    ref = seqs[0]

    # reference numbering map
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
            tri = s[col:col+3]

            if len(tri) < 3:
                continue
            if any(ch in IGNORE_AA for ch in tri):
                continue

            covered_n += 1
            if is_n_glyco_motif(tri):
                motif_n += 1

        coverage_frac = covered_n / n_seq if n_seq > 0 else 0.0

        if coverage_frac < min_triplet_coverage:
            continue

        motif_freq = motif_n / covered_n if covered_n > 0 else 0.0
        motif_pct = motif_freq * 100

        rows.append({
            "Position": start_pos,
            "Covered_n": covered_n,
            "Coverage_pct": round(coverage_frac * 100, 2),
            "Motif_n": motif_n,
            "Motif_freq": round(motif_freq, 4),
            "Motif_pct": round(motif_pct, 2)
        })

    return pd.DataFrame(rows)

def analyze_period(label, fasta_file, min_triplet_coverage, presence_threshold):
    names, seqs = read_aligned_fasta(fasta_file)
    df = compute_glyco_prevalence(seqs, min_triplet_coverage=min_triplet_coverage)

    safe_label = label.replace("-", "_")
    df.to_csv(f"{safe_label}_glyco_site_prevalence.csv", index=False)

    present_df = df[df["Motif_pct"] >= presence_threshold].copy()
    present_df.to_csv(f"{safe_label}_glyco_sites_present_{int(presence_threshold)}pct.csv", index=False)

    summary = {
        "Period": label,
        "Input_file": fasta_file,
        "Sequence_count": len(seqs),
        "Alignment_length": len(seqs[0]),
        "Analysed_positions": len(df),
        "Predicted_sites_ge_threshold": len(present_df),
        "Mean_motif_prevalence_pct": round(df["Motif_pct"].mean(), 2),
        "Max_motif_prevalence_pct": round(df["Motif_pct"].max(), 2) if len(df) > 0 else 0.0
    }

    return summary, df

def classify_trend(vals, threshold):
    """
    vals = [period1_pct, period2_pct, period3_pct]
    threshold = presence threshold in %
    """
    p = [v >= threshold for v in vals]

    if all(p):
        return "Persistent"
    elif (not p[0]) and p[-1]:
        return "Gain"
    elif p[0] and (not p[-1]):
        return "Loss"
    else:
        return "Transient"

def build_combined_matrix(period_labels, period_dfs, presence_threshold, min_display_threshold):
    all_positions = sorted(set().union(*[set(df["Position"]) for df in period_dfs.values()]))

    combined = pd.DataFrame({"Position": all_positions})

    for label in period_labels:
        temp = period_dfs[label][["Position", "Motif_pct", "Motif_n", "Covered_n", "Coverage_pct"]].copy()
        safe = label.replace("-", "_")
        temp.columns = [
            "Position",
            f"{safe}_pct",
            f"{safe}_motif_n",
            f"{safe}_covered_n",
            f"{safe}_coverage_pct"
        ]
        combined = combined.merge(temp, on="Position", how="left")

    combined = combined.fillna(0)

    pct_cols = [f"{label.replace('-', '_')}_pct" for label in period_labels]

    combined["Max_pct"] = combined[pct_cols].max(axis=1)
    combined["Min_pct"] = combined[pct_cols].min(axis=1)
    combined["Dynamic_range"] = combined["Max_pct"] - combined["Min_pct"]

    combined["Trend"] = combined.apply(
        lambda row: classify_trend([row[pct_cols[0]], row[pct_cols[1]], row[pct_cols[2]]], presence_threshold),
        axis=1
    )

    combined["Display"] = combined["Max_pct"] >= min_display_threshold

    # Sort for figure
    trend_order = {"Persistent": 0, "Gain": 1, "Loss": 2, "Transient": 3}
    combined["Trend_order"] = combined["Trend"].map(trend_order)
    combined = combined.sort_values(["Display", "Trend_order", "Max_pct", "Dynamic_range"], ascending=[False, True, False, False])

    displayed = combined[combined["Display"]].copy()

    if displayed.empty:
        print("[WARNING] No glycosylation sites reached the display threshold.")
        print("[WARNING] Lowering display threshold to 1%.")
        combined["Display"] = combined["Max_pct"] >= 1.0
        displayed = combined[combined["Display"]].copy()

    combined.to_csv("glyco_evolution_all_positions.csv", index=False)
    displayed.to_csv("glyco_evolution_displayed_positions.csv", index=False)

    return combined, displayed

def plot_glyco_evolution(displayed_df, summary_df, period_labels, out_prefix, title, presence_threshold, min_display_threshold, top_dynamic_n):
    safe_labels = [p.replace("-", "_") for p in period_labels]
    pct_cols = [f"{s}_pct" for s in safe_labels]

    heatmap_data = displayed_df[pct_cols].values
    site_labels = displayed_df["Position"].astype(int).astype(str).tolist()

    # top dynamic sites
    top_dynamic = displayed_df.sort_values(["Dynamic_range", "Max_pct"], ascending=False).head(top_dynamic_n).copy()

    trend_counts = displayed_df["Trend"].value_counts().reindex(["Persistent", "Gain", "Loss", "Transient"], fill_value=0)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10
    })

    fig = plt.figure(figsize=(16, 11), facecolor="white")

    # ---------------- Panel A: Heatmap ----------------
    ax1 = fig.add_axes([0.08, 0.48, 0.58, 0.42])
    cax = fig.add_axes([0.68, 0.48, 0.015, 0.42])

    im = ax1.imshow(heatmap_data, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)

    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Predicted glycosylation prevalence (%)")

    ax1.set_xticks(np.arange(len(period_labels)))
    ax1.set_xticklabels([f"{p}\n(n={n})" for p, n in zip(summary_df["Period"], summary_df["Sequence_count"])])
    ax1.set_yticks(np.arange(len(site_labels)))
    ax1.set_yticklabels(site_labels)
    ax1.set_xlabel("Sequence collection period")
    ax1.set_ylabel("NA protein position (motif start)")
    ax1.set_title(title, pad=14)
    ax1.text(-0.08, 1.03, "A", transform=ax1.transAxes, fontsize=15, fontweight="bold")

    # annotate heatmap values
    if len(site_labels) <= 18:
        for i in range(heatmap_data.shape[0]):
            for j in range(heatmap_data.shape[1]):
                val = heatmap_data[i, j]
                ax1.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=9, color="black")

    # ---------------- Panel B: Trend counts ----------------
    ax2 = fig.add_axes([0.74, 0.58, 0.22, 0.25])

    trend_colors = {
        "Persistent": "#2F6DB3",
        "Gain": "#2E8B57",
        "Loss": "#C62828",
        "Transient": "#BFBFBF"
    }

    cats = trend_counts.index.tolist()
    vals = trend_counts.values

    ax2.bar(cats, vals, color=[trend_colors[c] for c in cats], edgecolor="black", linewidth=0.6)

    for i, v in enumerate(vals):
        ax2.text(i, v + 0.15, str(v), ha="center", va="bottom", fontsize=10)

    ax2.set_ylabel("Number of sites")
    ax2.set_title("Trend category counts")
    ax2.grid(axis="y", linestyle="--", alpha=0.25)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(-0.20, 1.05, "B", transform=ax2.transAxes, fontsize=15, fontweight="bold")

    # ---------------- Panel C: Top dynamic site trajectories ----------------
    ax3 = fig.add_axes([0.08, 0.12, 0.88, 0.24])

    x = np.arange(len(period_labels))

    for _, row in top_dynamic.iterrows():
        pos = int(row["Position"])
        y = [row[pct_cols[0]], row[pct_cols[1]], row[pct_cols[2]]]
        ax3.plot(x, y, marker="o", linewidth=2, label=f"Site {pos}")

    ax3.set_xticks(x)
    ax3.set_xticklabels([f"{p}\n(n={n})" for p, n in zip(summary_df["Period"], summary_df["Sequence_count"])])
    ax3.set_ylabel("Predicted glycosylation prevalence (%)")
    ax3.set_xlabel("Sequence collection period")
    ax3.set_title(f"Top {min(top_dynamic_n, len(top_dynamic))} dynamic glycosylation sites")
    ax3.grid(axis="y", linestyle="--", alpha=0.25)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.legend(frameon=False, ncol=4, loc="upper left")
    ax3.text(-0.03, 1.05, "C", transform=ax3.transAxes, fontsize=15, fontweight="bold")

    # ---------------- Global legend / note ----------------
    legend_handles = [
        Patch(facecolor="#2F6DB3", edgecolor="black", label="Persistent site"),
        Patch(facecolor="#2E8B57", edgecolor="black", label="Gain site"),
        Patch(facecolor="#C62828", edgecolor="black", label="Loss site"),
        Patch(facecolor="#BFBFBF", edgecolor="black", label="Transient site")
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.72, 0.93),
        ncol=2,
        frameon=False
    )

    fig.text(
        0.08,
        0.04,
        f"Potential N-linked glycosylation motifs were defined as N-X-S/T (X ≠ P). "
        f"Site prevalence indicates the percentage of aligned sequences containing the motif at a given NA position. "
        f"Trend classification used a {presence_threshold:.1f}% prevalence threshold. "
        f"Only sites with maximum prevalence ≥ {min_display_threshold:.1f}% across at least one period are shown in the heatmap. "
        f"Positions with <80% triplet coverage were excluded. The 2021–2025 group has a smaller sample size and should be interpreted cautiously.",
        fontsize=9,
        style="italic"
    )

    png_file = f"{out_prefix}.png"
    pdf_file = f"{out_prefix}.pdf"
    svg_file = f"{out_prefix}.svg"

    plt.savefig(png_file, dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(pdf_file, bbox_inches="tight", facecolor="white")
    plt.savefig(svg_file, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"[DONE] Saved figure: {png_file}")
    print(f"[DONE] Saved figure: {pdf_file}")
    print(f"[DONE] Saved figure: {svg_file}")

def main():
    parser = argparse.ArgumentParser(description="Generate a journal-level glycosylation gain/loss evolution map for Asian H1N1 NA protein.")
    parser.add_argument("--f1", required=True, help="Aligned FASTA for 2015-2017")
    parser.add_argument("--f2", required=True, help="Aligned FASTA for 2018-2020")
    parser.add_argument("--f3", required=True, help="Aligned FASTA for 2021-2025")
    parser.add_argument("--l1", default="2015-2017")
    parser.add_argument("--l2", default="2018-2020")
    parser.add_argument("--l3", default="2021-2025")
    parser.add_argument("-o", "--out", default="H1N1_NA_469_glyco_evolution_map", help="Output prefix")
    parser.add_argument("--title", default="Evolution map of potential N-linked glycosylation motifs in Asian H1N1 neuraminidase (NA)", help="Figure title")
    parser.add_argument("--presence_threshold", type=float, default=5.0, help="Prevalence threshold (%) to classify gain/loss/persistent sites")
    parser.add_argument("--min_display_threshold", type=float, default=5.0, help="Minimum max prevalence (%) required for display in heatmap")
    parser.add_argument("--min_triplet_coverage", type=float, default=0.80, help="Minimum fraction of sequences with non-gap triplet coverage")
    parser.add_argument("--top_dynamic", type=int, default=8, help="Number of top dynamic sites to show in Panel C")
    args = parser.parse_args()

    period_labels = [args.l1, args.l2, args.l3]
    period_files = [args.f1, args.f2, args.f3]

    summaries = []
    period_dfs = {}

    for label, fasta in zip(period_labels, period_files):
        print(f"\n[INFO] Analysing {label}: {fasta}")
        summary, df = analyze_period(
            label=label,
            fasta_file=fasta,
            min_triplet_coverage=args.min_triplet_coverage,
            presence_threshold=args.presence_threshold
        )
        summaries.append(summary)
        period_dfs[label] = df

    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(f"{args.out}_period_summary.csv", index=False)

    combined_df, displayed_df = build_combined_matrix(
        period_labels=period_labels,
        period_dfs=period_dfs,
        presence_threshold=args.presence_threshold,
        min_display_threshold=args.min_display_threshold
    )

    top_dynamic_df = displayed_df.sort_values(["Dynamic_range", "Max_pct"], ascending=False).head(args.top_dynamic)
    top_dynamic_df.to_csv(f"{args.out}_top_dynamic_sites.csv", index=False)

    trend_summary = displayed_df["Trend"].value_counts().reindex(["Persistent", "Gain", "Loss", "Transient"], fill_value=0)
    trend_summary = trend_summary.reset_index()
    trend_summary.columns = ["Trend", "Count"]
    trend_summary.to_csv(f"{args.out}_trend_summary.csv", index=False)

    print("\n[PERIOD SUMMARY]")
    print(summary_df.to_string(index=False))

    print("\n[DISPLAYED GLYCOSYLATION SITES]")
    print(displayed_df[["Position", "Trend", "Max_pct", "Dynamic_range"]].head(20).to_string(index=False))

    plot_glyco_evolution(
        displayed_df=displayed_df,
        summary_df=summary_df,
        period_labels=period_labels,
        out_prefix=args.out,
        title=args.title,
        presence_threshold=args.presence_threshold,
        min_display_threshold=args.min_display_threshold,
        top_dynamic_n=args.top_dynamic
    )

if __name__ == "__main__":
    main()
