import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =========================
# INPUT / OUTPUT PATH
# =========================
base_dir = Path.home() / "Desktop" / "influenza_2"
csv_path = base_dir / "H1N1_NA_mutation_frequency.csv"
out_dir = base_dir / "report_figures"
out_dir.mkdir(exist_ok=True)

# =========================
# LOAD CSV
# =========================
df = pd.read_csv(csv_path)

# =========================
# AUTO-DETECT COLUMN NAMES
# =========================
possible_pos_cols = ["Position", "position", "NA_Position", "Site", "Residue"]
possible_count_cols = ["Mutation_Count", "mutation_count", "Count", "Frequency", "MutationFrequency"]

pos_col = None
count_col = None

for c in possible_pos_cols:
    if c in df.columns:
        pos_col = c
        break

for c in possible_count_cols:
    if c in df.columns:
        count_col = c
        break

if pos_col is None or count_col is None:
    raise ValueError(
        f"Column detect korte pari nai.\nAvailable columns: {list(df.columns)}\n"
        "Ensure there is one position column and one mutation count column."
    )

df = df[[pos_col, count_col]].copy()
df.columns = ["Position", "Mutation_Count"]

# numeric kore nao
df["Position"] = pd.to_numeric(df["Position"], errors="coerce")
df["Mutation_Count"] = pd.to_numeric(df["Mutation_Count"], errors="coerce")
df = df.dropna().sort_values("Position").reset_index(drop=True)

# =========================
# TOP PEAKS SELECT
# =========================
top_n = 12
top_df = df.nlargest(top_n, "Mutation_Count").copy()

# majhkhane important peak miss hole force include
force_positions = [188, 314]
forced_df = df[df["Position"].isin(force_positions)].copy()

top_df = pd.concat([top_df, forced_df], ignore_index=True)
top_df = top_df.drop_duplicates(subset=["Position"])
top_df = top_df.sort_values("Position").reset_index(drop=True)

# save labeled peak list
top_df.to_csv(out_dir / "Top_peak_residues_labeled_final.csv", index=False)

# =========================
# FIGURE
# =========================
plt.rcParams["font.family"] = "DejaVu Sans"

fig, ax = plt.subplots(figsize=(16, 8))

# landscape line plot
ax.plot(
    df["Position"],
    df["Mutation_Count"],
    color="#2F6DA3",
    linewidth=1.8
)

# top peaks mark
ax.scatter(
    top_df["Position"],
    top_df["Mutation_Count"],
    color="#2F6DA3",
    s=14,
    zorder=4
)

# =========================
# LABEL OFFSETS
# close to peaks, but not overlapping
# =========================
offset_patterns = [
    (-8, 8),
    (0, 8),
    (8, 8),
    (-10, 18),
    (10, 18),
    (-12, 28),
    (12, 28),
    (-6, 16),
    (6, 16),
]

last_x = -999
cluster_index = 0

for _, row in top_df.iterrows():
    x = int(row["Position"])
    y = float(row["Mutation_Count"])

    # jodi पासাপासी peak hoy, different offset use korbe
    if x - last_x <= 12:
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
        fontsize=11,
        fontweight="bold",
        color="black",
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.82,
            pad=0.5
        )
    )

    last_x = x

# =========================
# STYLE
# =========================
ax.set_title(
    "Mutation frequency landscape of H1N1 Neuraminidase gene (Asia, 2015-2025)",
    fontsize=18,
    fontweight="bold",
    pad=14
)

ax.set_xlabel("NA position", fontsize=14)
ax.set_ylabel("Mutation count", fontsize=14)

ax.tick_params(axis="x", labelsize=10)
ax.tick_params(axis="y", labelsize=10)

ax.grid(axis="y", linestyle="--", alpha=0.25)

# x range একটু buffer
xmin = max(0, df["Position"].min() - 5)
xmax = df["Position"].max() + 5
ax.set_xlim(xmin, xmax)

# y buffer
ymax = df["Mutation_Count"].max()
ax.set_ylim(0, ymax * 1.18)

plt.tight_layout()

# =========================
# SAVE
# =========================
png_out = out_dir / "Figure_mutation_frequency_landscape_final_clean.png"
pdf_out = out_dir / "Figure_mutation_frequency_landscape_final_clean.pdf"
svg_out = out_dir / "Figure_mutation_frequency_landscape_final_clean.svg"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")
plt.savefig(svg_out, bbox_inches="tight", facecolor="white")

print("[DONE] Figure saved:")
print(png_out)
print(pdf_out)
print(svg_out)

print("\n[DONE] Top peak residue list saved:")
print(out_dir / "Top_peak_residues_labeled_final.csv")
