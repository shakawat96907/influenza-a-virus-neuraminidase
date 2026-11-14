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

# column clean
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
# TOP PEAKS
# =========================
top_n = 12
top_df = df.nlargest(top_n, "Mutation_Count").sort_values("Position").reset_index(drop=True)
top_df.to_csv(out_dir / "Top_peak_residues_clean_labels_close.csv", index=False)

# =========================
# FIGURE STYLE
# =========================
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.linewidth": 1.0
})

fig, ax = plt.subplots(figsize=(11.69, 8.27))   # A4 landscape

# main landscape line
ax.plot(
    positions,
    counts,
    color="#2E6FBB",
    linewidth=1.6
)

# peak dots
ax.scatter(
    top_df["Position"],
    top_df["Mutation_Count"],
    color="#1F4E79",
    s=18,
    zorder=3
)

# title and labels
ax.set_title(
    "Mutation frequency landscape of H1N1 Neuraminidase gene (Asia, 2015-2025)",
    fontsize=16,
    fontweight="bold",
    pad=14
)

ax.set_xlabel("NA position", fontsize=13)
ax.set_ylabel("Mutation count", fontsize=13)

ax.set_xlim(1, max(positions) + 5)
ax.set_ylim(0, max(counts) * 1.15)

ax.set_xticks(np.arange(0, max(positions) + 1, 25))
ax.tick_params(axis="x", labelsize=9)
ax.tick_params(axis="y", labelsize=10)

ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.35)

# =========================
# CLOSE BUT CLEAN LABELS
# =========================
# peak er kache label thakbe
# clustered peak হলে ছোট ছোট offset change হবে
offset_patterns = [
    (0, 8),     # just above
    (6, 10),    # slightly right
    (-6, 10),   # slightly left
    (0, 14),    # একটু ওপরে
    (8, 14),    # right-up
    (-8, 14),   # left-up
    (0, 18)     # more up if needed
]

last_x = -999
cluster_index = 0

for _, row in top_df.iterrows():
    x = int(row["Position"])
    y = float(row["Mutation_Count"])

    # যদি আগের peak এর খুব কাছে হয়, অন্য offset use করবে
    if x - last_x <= 10:
        cluster_index = (cluster_index + 1) % len(offset_patterns)
    else:
        cluster_index = 0

    dx, dy = offset_patterns[cluster_index]

    ax.annotate(
        str(x),
        xy=(x, y),
        xytext=(dx, dy),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color="black",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.80,
            pad=0.6
        )
    )

    last_x = x

plt.tight_layout()

# =========================
# SAVE
# =========================
png_out = out_dir / "Figure_mutation_frequency_landscape_clean_labels_close.png"
pdf_out = out_dir / "Figure_mutation_frequency_landscape_clean_labels_close.pdf"
svg_out = out_dir / "Figure_mutation_frequency_landscape_clean_labels_close.svg"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.savefig(svg_out, bbox_inches="tight", facecolor="white")

print("[DONE] Saved:")
print(png_out)
print(pdf_out)
print(svg_out)
print(out_dir / "Top_peak_residues_clean_labels_close.csv")
