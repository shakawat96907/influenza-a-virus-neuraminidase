from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

# Read aligned fasta
alignment = AlignIO.read("H1N1_NA_Asia_aligned.fasta", "fasta")

# Reference sequence
ref = alignment[0].seq
mut = {i: Counter() for i in range(len(ref))}

# Count mutations at each position
for rec in alignment[1:]:
    for i, aa in enumerate(rec.seq):
        if aa != ref[i] and aa != "-":
            mut[i][aa] += 1

# Build mutation table
rows = []
for pos, c in mut.items():
    for aa, count in c.items():
        rows.append([pos + 1, ref[pos], aa, count])

df = pd.DataFrame(rows, columns=["Position", "Reference", "Mutation", "Count"])
df.to_csv("H1N1_NA_mutation_frequency.csv", index=False)

print("Mutation frequency file created")

# Total mutation count per position
total_df = df.groupby("Position", as_index=False)["Count"].sum()

# Dominant mutation per position
dominant_idx = df.groupby("Position")["Count"].idxmax()
dominant_df = df.loc[dominant_idx, ["Position", "Reference", "Mutation", "Count"]].copy()

# Make dominant mutation label like N248D
dominant_df["Dominant_Label"] = (
    dominant_df["Reference"].astype(str)
    + dominant_df["Position"].astype(str)
    + dominant_df["Mutation"].astype(str)
)

# Merge total counts with dominant labels
plot_df = pd.merge(
    total_df,
    dominant_df[["Position", "Dominant_Label"]],
    on="Position",
    how="left"
)

# Select top 15 positions by total mutation count
plot_df = plot_df.sort_values("Count", ascending=False).head(15)

# Sort by count for plotting
plot_df = plot_df.sort_values("Count", ascending=False)

# Create publication-style bar plot
plt.figure(figsize=(12, 6))
bars = plt.bar(plot_df["Dominant_Label"], plot_df["Count"], edgecolor="black", linewidth=0.8)

plt.xlabel("Dominant Mutation", fontsize=12)
plt.ylabel("Mutation Count", fontsize=12)
plt.title("Top 15 Dominant H1N1 NA Mutations in Asia", fontsize=14)

plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(fontsize=10)

# Add count values on top of bars
for bar, count in zip(bars, plot_df["Count"]):
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 1,
        str(count),
        ha="center",
        va="bottom",
        fontsize=9
    )

plt.tight_layout()
plt.savefig("H1N1_NA_top15_publication_plot.png", dpi=300, bbox_inches="tight")
plt.show()

print("Top 15 publication-style bar plot created")
