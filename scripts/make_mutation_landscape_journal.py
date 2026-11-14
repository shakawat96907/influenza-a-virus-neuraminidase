import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# =========================
# Input / Output
# =========================
input_file = "H1N1_NA_mutation_frequency.csv"
out_dir = Path("report_figures")
out_dir.mkdir(exist_ok=True)

# =========================
# Read CSV
# =========================
df = pd.read_csv(input_file)

# column names clean
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
print("Detected columns:", df.columns.tolist())

# auto-detect possible column names
position_candidates = ["position", "site", "residue", "residue_position", "na_position", "pos"]
count_candidates = ["mutation_count", "count", "frequency", "mutation_frequency", "mut_count", "mutations"]

pos_col = next((c for c in position_candidates if c in df.columns), None)
count_col = next((c for c in count_candidates if c in df.columns), None)

if pos_col is None:
    raise ValueError("Position column not found in CSV.")
if count_col is None:
    raise ValueError("Mutation count column not found in CSV.")

df = df[[pos_col, count_col]].copy()
df = df.rename(columns={pos_col: "Position", count_col: "Mutation_Count"})

df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
df["Mutation_Count"] = pd.to_numeric(df["Mutation_Count"], errors="coerce")
df = df.dropna().sort_values("Position").reset_index(drop=True)

positions = df["Position"].astype(int).values
counts = df["Mutation_Count"].values

# =========================
# Select top peaks to label
# =========================
top_n = 15
top_df = df.nlargest(top_n, "Mutation_Count").sort_values("Position").reset_index(drop=True)
top_df.to_csv(out_dir / "Top_peak_residues_landscape.csv", index=False)

# =========================
# Plot style
# =========================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 1.0
})

fig, ax = plt.subplots(figsize=(11.69, 8.27))   # A4 landscape

# main line plot
ax.plot(
    positions,
    counts,
    color="#2E6FBB",
    linewidth=1.8
)

# small peak markers for top labeled residues
ax.scatter(
    top_df["Position"],
    top_df["Mutation_Count"],
    color="#1F4E79",
    s=18,
    zorder=3
)

# =========================
# Labels and title
# =========================
ax.set_title(
    "Mutation frequency landscape of H1N1 Neuraminidase gene (Asia, 2015-2025)",
    fontsize=16,
    fontweight="bold",
    pad=14
)

ax.set_xlabel("NA position", fontsize=13)
ax.set_ylabel("Mutation count", fontsize=13)

# =========================
# Axis formatting
# =========================
ax.set_xlim(1, max(positions) + 5)
ax.set_ylim(0, max(counts) * 1.22)

xticks = np.arange(0, max(positions) + 1, 25)
ax.set_xticks(xticks)
ax.tick_params(axis="x", labelsize=9)
ax.tick_params(axis="y", labelsize=10)

ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.35)
ax.spines["top"].set_visible(True)
ax.spines["right"].set_visible(True)

# =========================
# Clean peak annotations
# =========================
last_x = -999
level = 0
ymax = max(counts)

# staggered offsets so labels don't overlap
offset_levels = [0.04, 0.08, 0.12, 0.16, 0.20]
x_shift_levels = [0, -4, 4, -7, 7]

for _, row in top_df.iterrows():
    x = int(row["Position"])
    y = float(row["Mutation_Count"])

    if x - last_x < 14:
        level = (level + 1) % len(offset_levels)
    else:
        level = 0

    y_text = y + ymax * offset_levels[level]
    x_text = x + x_shift_levels[level]

    ax.annotate(
        str(x),
        xy=(x, y),
        xytext=(x_text, y_text),
        textcoords="data",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color="black",
        arrowprops=dict(
            arrowstyle="-",
            color="gray",
            lw=0.7,
            shrinkA=0,
            shrinkB=0
        ),
        bbox=dict(
            boxstyle="round,pad=0.15",
            fc="white",
            ec="none",
            alpha=0.85
        )
    )

    last_x = x

# =========================
# Save files
# =========================
plt.tight_layout()

png_out = out_dir / "Figure_mutation_frequency_landscape_H1N1_Asia_2015_2025.png"
pdf_out = out_dir / "Figure_mutation_frequency_landscape_H1N1_Asia_2015_2025.pdf"
svg_out = out_dir / "Figure_mutation_frequency_landscape_H1N1_Asia_2015_2025.svg"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.savefig(svg_out, bbox_inches="tight", facecolor="white")

print("[DONE] Saved:")
print(png_out)
print(pdf_out)
print(svg_out)
print(out_dir / "Top_peak_residues_landscape.csv")
