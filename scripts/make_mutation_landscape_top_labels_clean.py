import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# =========================
# Input / output
# =========================
input_file = "H1N1_NA_mutation_frequency.csv"
out_dir = Path("report_figures")
out_dir.mkdir(exist_ok=True)

TOP_N = 15   # koyta top peak label dibe

# =========================
# Read CSV
# =========================
df = pd.read_csv(input_file)
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

print("Detected columns:", df.columns.tolist())

# Auto-detect columns
position_candidates = ["position", "site", "residue", "residue_position", "na_position", "pos"]
count_candidates = ["mutation_count", "count", "frequency", "mutation_frequency", "mut_count", "mutations"]

pos_col = next((c for c in position_candidates if c in df.columns), None)
count_col = next((c for c in count_candidates if c in df.columns), None)

if pos_col is None:
    raise ValueError("Position column not found.")
if count_col is None:
    raise ValueError("Mutation count/frequency column not found.")

df = df[[pos_col, count_col]].copy()
df.columns = ["Position", "Mutation_Count"]

df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
df["Mutation_Count"] = pd.to_numeric(df["Mutation_Count"], errors="coerce")
df = df.dropna().sort_values("Position")

# =========================
# Select top peaks
# =========================
top_df = df.nlargest(TOP_N, "Mutation_Count").copy()
top_df = top_df.sort_values("Position").reset_index(drop=True)
top_df.to_csv(out_dir / "Top_peaks_labelled_positions.csv", index=False)

# =========================
# Plot style
# =========================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.9,
})

fig, ax = plt.subplots(figsize=(14, 7))  # landscape

# Line plot
ax.plot(
    df["Position"],
    df["Mutation_Count"],
    color="#3B73B9",
    linewidth=2.0
)

# Small dots only on labelled peaks
ax.scatter(
    top_df["Position"],
    top_df["Mutation_Count"],
    color="#1E4E8C",
    s=24,
    zorder=3
)

# =========================
# Title and axis
# =========================
ax.set_title(
    "Mutation Landscape of H1N1 Neuraminidase (Asia, 2015–2025)",
    fontsize=20,
    pad=14
)

ax.set_xlabel("NA Position", fontsize=16)
ax.set_ylabel("Mutation Count", fontsize=16)

ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)

ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.28)

# =========================
# Non-overlapping peak labels
# =========================
ymax = df["Mutation_Count"].max()
label_levels = [0.035, 0.085, 0.135, 0.185]  # staggered heights

last_x = -999
level_index = 0

for _, row in top_df.iterrows():
    x = int(row["Position"])
    y = row["Mutation_Count"]

    # If peaks are close, move label to next level
    if abs(x - last_x) < 12:
        level_index = (level_index + 1) % len(label_levels)
    else:
        level_index = 0

    y_text = y + ymax * label_levels[level_index]

    ax.text(
        x,
        y_text,
        str(x),       # only residue number
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        rotation=90,
        color="#222222",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.8)
    )

    # short connector line
    ax.plot(
        [x, x],
        [y + ymax * 0.01, y_text - ymax * 0.01],
        color="#555555",
        linewidth=0.5,
        alpha=0.55
    )

    last_x = x

# Axis limits
ax.set_xlim(0, df["Position"].max() + 8)
ax.set_ylim(0, ymax * 1.25)

plt.tight_layout()

# =========================
# Save
# =========================
png_out = out_dir / "H1N1_NA_mutation_landscape_top_labels_clean.png"
pdf_out = out_dir / "H1N1_NA_mutation_landscape_top_labels_clean.pdf"
svg_out = out_dir / "H1N1_NA_mutation_landscape_top_labels_clean.svg"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.savefig(svg_out, bbox_inches="tight", facecolor="white")

print("[DONE] Saved:")
print(png_out)
print(pdf_out)
print(svg_out)
print(out_dir / "Top_peaks_labelled_positions.csv")
