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
    raise ValueError("Mutation count column not found.")

df = df[[pos_col, count_col]].copy()
df = df.rename(columns={pos_col: "Position", count_col: "Mutation_Count"})

df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
df["Mutation_Count"] = pd.to_numeric(df["Mutation_Count"], errors="coerce")
df = df.dropna().sort_values("Position").reset_index(drop=True)

positions = df["Position"].astype(int).values
counts = df["Mutation_Count"].values

# =========================
# Select top peak residues
# =========================
top_n = 15
top_df = df.nlargest(top_n, "Mutation_Count").sort_values("Position").reset_index(drop=True)
top_df.to_csv(out_dir / "Top15_overall_mutation_frequency_sites.csv", index=False)

# =========================
# A4 landscape plot
# =========================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 0.9,
})

fig, ax = plt.subplots(figsize=(11.69, 8.27))  # A4 landscape

ax.bar(
    positions,
    counts,
    width=0.85,
    color="#4A90E2",
    edgecolor="#2F5D8A",
    linewidth=0.35
)

# =========================
# Title and axis labels
# =========================
ax.set_title(
    "H1N1 NA Mutation Frequencies (All Positions)",
    fontsize=17,
    fontweight="bold",
    pad=14
)

ax.set_xlabel("NA Position", fontsize=14)
ax.set_ylabel("Mutation Count", fontsize=14)

# =========================
# Axis formatting
# =========================
ax.set_xlim(0, max(positions) + 5)
ax.set_ylim(0, max(counts) * 1.25)

xticks = np.arange(0, max(positions) + 1, 25)
ax.set_xticks(xticks)

ax.tick_params(axis="x", labelsize=8, rotation=90)
ax.tick_params(axis="y", labelsize=11)

ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.35)
ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# =========================
# Non-overlapping top peak labels
# =========================
last_x = -999
level = 0
levels = [0.035, 0.075, 0.115, 0.155]  # staggered heights

for _, row in top_df.iterrows():
    x = int(row["Position"])
    y = row["Mutation_Count"]

    # nearby peak hole label height alternate korbe
    if x - last_x < 12:
        level = (level + 1) % len(levels)
    else:
        level = 0

    label_y = y + max(counts) * levels[level]

    ax.text(
        x,
        label_y,
        str(x),
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
        rotation=90,
        color="#1F1F1F",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.75,
            pad=0.8
        )
    )

    last_x = x

# =========================
# Save outputs
# =========================
plt.tight_layout()

png_out = out_dir / "Figure_overall_mutation_frequency_A4_landscape_no_overlap.png"
pdf_out = out_dir / "Figure_overall_mutation_frequency_A4_landscape_no_overlap.pdf"
svg_out = out_dir / "Figure_overall_mutation_frequency_A4_landscape_no_overlap.svg"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.savefig(svg_out, bbox_inches="tight", facecolor="white")

print("[DONE] Saved:")
print(png_out)
print(pdf_out)
print(svg_out)
print(out_dir / "Top15_overall_mutation_frequency_sites.csv")
