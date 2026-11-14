#!/usr/bin/env python3

import argparse
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from Bio import SeqIO


IGNORE_AA = set(["-", "X", "?", ".", "*"])


def read_aligned_protein_fasta(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if len(records) == 0:
        raise ValueError("No FASTA records found.")

    names = []
    seqs = []

    for rec in records:
        seq = str(rec.seq).strip().upper().replace(" ", "")
        names.append(rec.id)
        seqs.append(seq)

    lengths = set(len(s) for s in seqs)
    if len(lengths) != 1:
        raise ValueError(
            "Sequences are not the same length. Please use an ALIGNED protein FASTA file."
        )

    return names, seqs


def compute_site_metrics(seqs, conserved_cutoff=0.02, hyper_cutoff=0.10, min_non_gap_frac=0.80):
    n_seq = len(seqs)
    aln_len = len(seqs[0])

    ref = seqs[0]
    rows = []
    ref_pos = 0

    for col in range(aln_len):
        ref_aa = ref[col]

        # reference numbering follows first sequence non-gap positions
        if ref_aa in IGNORE_AA:
            continue
        ref_pos += 1

        col_aas = [s[col] for s in seqs if s[col] not in IGNORE_AA]
        non_gap_count = len(col_aas)
        non_gap_frac = non_gap_count / n_seq

        if non_gap_frac < min_non_gap_frac:
            continue

        counts = Counter(col_aas)
        consensus_aa, consensus_count = counts.most_common(1)[0]

        conservation = consensus_count / non_gap_count
        variation = 1.0 - conservation

        # Shannon entropy
        probs = np.array([c / non_gap_count for c in counts.values()])
        entropy = -np.sum(probs * np.log2(probs))

        # top alternative residue relative to reference residue
        alt_counts = {aa: c for aa, c in counts.items() if aa != ref_aa}
        if alt_counts:
            top_alt, top_alt_count = sorted(alt_counts.items(), key=lambda x: x[1], reverse=True)[0]
            top_alt_freq = top_alt_count / non_gap_count
            mutation_label = f"{ref_aa}{ref_pos}{top_alt}"
        else:
            top_alt, top_alt_count, top_alt_freq, mutation_label = "", 0, 0.0, ""

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
            "Top_alt_frequency": round(top_alt_freq, 4),
            "Mutation_label": mutation_label,
            "Category": category
        })

    df = pd.DataFrame(rows)
    return df


def make_region_blocks(df):
    if df.empty:
        return pd.DataFrame(columns=["Start", "End", "Category", "Length"])

    blocks = []

    current_cat = df.iloc[0]["Category"]
    start = int(df.iloc[0]["Position"])
    prev = start

    for _, row in df.iloc[1:].iterrows():
        pos = int(row["Position"])
        cat = row["Category"]

        if cat == current_cat and pos == prev + 1:
            prev = pos
        else:
            blocks.append({
                "Start": start,
                "End": prev,
                "Category": current_cat,
                "Length": prev - start + 1
            })
            start = pos
            prev = pos
            current_cat = cat

    blocks.append({
        "Start": start,
        "End": prev,
        "Category": current_cat,
        "Length": prev - start + 1
    })

    return pd.DataFrame(blocks)


def plot_linear_landscape(
    df,
    blocks,
    out_prefix,
    title,
    conserved_cutoff,
    hyper_cutoff,
    top_n=12
):
    colors = {
        "Conserved": "#2F6DB3",      # blue
        "Variable": "#BFBFBF",       # grey
        "Hypervariable": "#C62828"   # red
    }

    max_pos = int(df["Position"].max())

    # sort top sites
    top_df = df.sort_values(["Variation", "Top_alt_frequency"], ascending=False).head(top_n).copy()
    hyper_df = df[df["Category"] == "Hypervariable"].copy()

    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "legend.fontsize": 10
    })

    fig = plt.figure(figsize=(16, 8), facecolor="white")
    gs = GridSpec(2, 1, height_ratios=[3.8, 1.2], hspace=0.18)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # -------- TOP PANEL: variation score --------
    x = df["Position"].values
    y = df["Variation"].values

    ax1.fill_between(x, 0, y, color="#D9D9D9", alpha=0.65, linewidth=0)
    ax1.plot(x, y, color="#4D4D4D", linewidth=1.3, label="Per-site variation score")

    if not hyper_df.empty:
        ax1.scatter(
            hyper_df["Position"],
            hyper_df["Variation"],
            color=colors["Hypervariable"],
            s=18,
            zorder=4,
            label="Hypervariable sites"
        )

    ax1.axhline(conserved_cutoff, color=colors["Conserved"], linestyle="--", linewidth=1.1)
    ax1.axhline(hyper_cutoff, color=colors["Hypervariable"], linestyle="--", linewidth=1.1)

    # annotate top sites
    if not top_df.empty:
        y_levels = [0.72, 0.80, 0.88, 0.96]
        for i, (_, row) in enumerate(top_df.iterrows()):
            pos = int(row["Position"])
            var = float(row["Variation"])
            label = row["Mutation_label"] if row["Mutation_label"] else f"Pos{pos}"
            text_y = y_levels[i % len(y_levels)]

            ax1.plot([pos, pos], [var, text_y - 0.015], color="black", linewidth=0.7, alpha=0.7)
            ax1.scatter([pos], [var], color="black", s=9, zorder=5)
            ax1.text(
                pos,
                text_y,
                label,
                ha="center",
                va="bottom",
                fontsize=9,
                rotation=45
            )

    ax1.set_xlim(1, max_pos)
    ax1.set_ylim(0, 1.02)
    ax1.set_ylabel("Variation score\n(1 - consensus frequency)")
    ax1.set_title(title, pad=14)
    ax1.grid(axis="y", linestyle="--", alpha=0.25)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # panel label
    ax1.text(-0.04, 1.03, "A", transform=ax1.transAxes, fontsize=14, fontweight="bold")

    # -------- BOTTOM PANEL: region bar --------
    y_bar = 0.25
    bar_height = 0.45

    for _, row in blocks.iterrows():
        start = int(row["Start"])
        end = int(row["End"])
        cat = row["Category"]

        ax2.broken_barh(
            [(start, end - start + 1)],
            (y_bar, bar_height),
            facecolors=colors[cat],
            edgecolors="black",
            linewidth=0.25
        )

    # draw protein baseline
    ax2.hlines(y_bar + bar_height / 2, 1, max_pos, color="black", linewidth=0.5, alpha=0.5)

    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.set_xlabel("NA protein amino acid position")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.grid(axis="x", linestyle="--", alpha=0.15)

    ax2.text(-0.04, 1.05, "B", transform=ax2.transAxes, fontsize=14, fontweight="bold")

    # -------- legend --------
    legend_handles = [
        Patch(facecolor=colors["Conserved"], edgecolor="black", label="Conserved region"),
        Patch(facecolor=colors["Variable"], edgecolor="black", label="Variable region"),
        Patch(facecolor=colors["Hypervariable"], edgecolor="black", label="Hypervariable region"),
        Line2D([0], [0], color="#4D4D4D", lw=1.3, label="Per-site variation score"),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors["Hypervariable"],
               markersize=6, label="Hypervariable site")
    ]

    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.99),
        ncol=5,
        frameon=False
    )

    plt.tight_layout(rect=[0.02, 0.03, 0.98, 0.95])

    png_file = f"{out_prefix}_linear_landscape.png"
    pdf_file = f"{out_prefix}_linear_landscape.pdf"
    svg_file = f"{out_prefix}_linear_landscape.svg"

    plt.savefig(png_file, dpi=600, bbox_inches="tight", facecolor="white")
    plt.savefig(pdf_file, bbox_inches="tight", facecolor="white")
    plt.savefig(svg_file, bbox_inches="tight", facecolor="white")
    plt.close()

    print(f"[DONE] Saved figure: {png_file}")
    print(f"[DONE] Saved figure: {pdf_file}")
    print(f"[DONE] Saved figure: {svg_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate publication-quality conserved-variable-hypervariable linear landscape plot from aligned H1N1 NA protein FASTA."
    )

    parser.add_argument("-i", "--input", required=True, help="Aligned protein FASTA")
    parser.add_argument("-o", "--out", default="H1N1_NA_Asia_469", help="Output prefix")
    parser.add_argument("--title", default="Linear variability landscape of Asian H1N1 neuraminidase (NA) protein", help="Figure title")
    parser.add_argument("--conserved_cutoff", type=float, default=0.02, help="Variation <= this value = conserved")
    parser.add_argument("--hyper_cutoff", type=float, default=0.10, help="Variation >= this value = hypervariable")
    parser.add_argument("--min_non_gap_frac", type=float, default=0.80, help="Minimum non-gap fraction per site")
    parser.add_argument("--top", type=int, default=12, help="Number of top variable sites to label")

    args = parser.parse_args()

    names, seqs = read_aligned_protein_fasta(args.input)

    print(f"[INFO] Number of sequences: {len(seqs)}")
    print(f"[INFO] Alignment length: {len(seqs[0])}")
    print(f"[INFO] Reference sequence for numbering: {names[0]}")

    df = compute_site_metrics(
        seqs,
        conserved_cutoff=args.conserved_cutoff,
        hyper_cutoff=args.hyper_cutoff,
        min_non_gap_frac=args.min_non_gap_frac
    )

    if df.empty:
        raise ValueError("No usable positions after filtering.")

    blocks = make_region_blocks(df)

    site_csv = f"{args.out}_site_variability.csv"
    blocks_csv = f"{args.out}_region_blocks.csv"
    top_csv = f"{args.out}_top_variable_sites.csv"

    df.to_csv(site_csv, index=False)
    blocks.to_csv(blocks_csv, index=False)
    df.sort_values(["Variation", "Top_alt_frequency"], ascending=False).head(args.top).to_csv(top_csv, index=False)

    print(f"[DONE] Saved site table: {site_csv}")
    print(f"[DONE] Saved region table: {blocks_csv}")
    print(f"[DONE] Saved top-site table: {top_csv}")

    print("\n[SUMMARY: CATEGORY COUNTS]")
    print(df["Category"].value_counts())

    print("\n[TOP VARIABLE SITES]")
    print(df.sort_values(["Variation", "Top_alt_frequency"], ascending=False).head(args.top)[[
        "Position", "Reference_AA", "Consensus_AA", "Variation",
        "Entropy", "Top_alt_AA", "Top_alt_frequency",
        "Mutation_label", "Category"
    ]].to_string(index=False))

    plot_linear_landscape(
        df=df,
        blocks=blocks,
        out_prefix=args.out,
        title=args.title,
        conserved_cutoff=args.conserved_cutoff,
        hyper_cutoff=args.hyper_cutoff,
        top_n=args.top
    )


if __name__ == "__main__":
    main()
