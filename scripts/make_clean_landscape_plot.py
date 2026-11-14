import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# =========================
# INPUT / OUTPUT
# =========================
input_file = "H1N1_NA_mutation_frequency.csv"
out_dir = Path("report_figures")
out_dir.mkdir(exist_ok=True)

# =========================
# READ CSV
# =========================
df = pd.read_csv(input_file)

# column name clean
df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]
print("Detected columns:", df.columns.tolist())

# possible column names
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
# SELECT MAJOR PEAKS
# =========================
top_n = 12   # চাইলে 10 বা 12 রাখতে পারো
top_df = df.nlargest(top_n, "Mutation_Count").sort_values("Position").reset_index(drop=True)
top_df.to_csv(out_dir / "Top_peak_residues_clean_labels.csv", index=False)

# =========================
# PLOT
# =========================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 1.0
})

fig, ax = plt.subplots(figsize=(11.69, 8.27))  # A4 landscape

# main line plot
ax.plot(
    positions,
    counts,
    color="#2E6FBB",
    linewidth=1.8
)

# peak markers
ax.scatter(
    top_df["Position"],
    top_df["Mutation_Count"],
    color="#1F4E79",
    s=20,
    zorder=3
)

# title and axis
ax.set_title(
    "Mutation frequency landscape of H1N1 Neuraminidase gene (Asia, 2015-2025)",
    fontsize=16,
    fontweight="bold",
    pad=14
)

ax.set_xlabel("NA position", fontsize=13)
ax.set_ylabel("Mutation count", fontsize=13)

ax.set_xlim(1, max(positions) + 5)
ax.set_ylim(0, max(counts) * 1.18)

ax.set_xticks(np.arange(0, max(positions) + 1, 25))
ax.tick_params(axis="x", labelsize=9)
ax.tick_params(axis="y", labelsize=10)

ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.35)

# =========================
# CLEAN LABELS (NO ARROWS / NO DAG)
# =========================
ymax = max(counts)

# alternate vertical offsets
y_offsets = [0.03, 0.07, 0.11, 0.15]
x_offsets = [-6, 6, -10, 10]

last_x = -999
level = 0

for _, row in top_df.iterrows():
    x = int(row["Position"])
    y = float(row["Mutation_Count"])

    if x - last_x < 18:
        level = (level + 1) % len(y_offsets)
    else:
        level = 0

    label_x = x + x_offsets[level]
    label_y = y + ymax * y_offsets[level]

    ax.text(
        label_x,
        label_y,
        str(x),
        fontsize=10,
        fontweight="bold",
        ha="center",
        va="bottom",
        color="black",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.85,
            pad=1.5
        )
    )

    last_x = x

# margins
plt.tight_layout()

# save
png_out = out_dir / "Figure_mutation_frequency_landscape_clean_labels.png"
pdf_out = out_dir / "Figure_mutation_frequency_landscape_clean_labels.pdf"
svg_out = out_dir / "Figure_mutation_frequency_landscape_clean_labels.svg"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.savefig(svg_out, bbox_inches="tight", facecolor="white")

print("[DONE] Saved:")
print(png_out)
print(pdf_out)
print(svg_out)
print(out_dir / "Top_peak_residues_clean_labels.csv")
