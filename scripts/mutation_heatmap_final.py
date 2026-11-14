import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load mutation frequency file
df = pd.read_csv("H1N1_NA_mutation_frequency.csv")

# Create mutation label (e.g., N248D)
df["Mutation_Label"] = (
    df["Reference"].astype(str)
    + df["Position"].astype(str)
    + df["Mutation"].astype(str)
)

# Select top 20 dominant mutations
top_n = 20
top_df = df.sort_values("Count", ascending=False).head(top_n).copy()

# Sort by mutation frequency
top_df = top_df.sort_values(["Count", "Position"], ascending=[False, True])

# Create heatmap matrix
heatmap_data = top_df.pivot_table(
    index="Mutation_Label",
    columns="Position",
    values="Count",
    fill_value=0
)

# Order rows by mutation frequency
row_order = top_df["Mutation_Label"].drop_duplicates().tolist()
heatmap_data = heatmap_data.reindex(row_order)

# Order columns by NA position
heatmap_data = heatmap_data.reindex(sorted(heatmap_data.columns), axis=1)

data = heatmap_data.values

# Plot figure
fig, ax = plt.subplots(figsize=(9,6))

im = ax.imshow(
    data,
    aspect="auto",
    cmap="YlOrRd",
    interpolation="nearest"
)

# Colorbar
cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Mutation Count", fontsize=11)

# Title (your requested title)
ax.set_title(
    "Heatmap of Top 20 Dominant Amino Acid Mutations in H1N1 Neuraminidase (Asia)",
    fontsize=14
)

# Axis labels
ax.set_xlabel("NA Position", fontsize=12)
ax.set_ylabel("Dominant Mutation", fontsize=12)

# X ticks
ax.set_xticks(np.arange(len(heatmap_data.columns)))
ax.set_xticklabels(
    heatmap_data.columns,
    rotation=45,
    ha="right",
    fontsize=8
)

# Y ticks
ax.set_yticks(np.arange(len(heatmap_data.index)))
ax.set_yticklabels(
    heatmap_data.index,
    fontsize=9
)

# Add cell borders
ax.set_xticks(np.arange(-0.5, len(heatmap_data.columns), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(heatmap_data.index), 1), minor=True)
ax.grid(which="minor", color="white", linestyle="-", linewidth=0.7)
ax.tick_params(which="minor", bottom=False, left=False)

# Add mutation count annotation
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        value = data[i, j]
        if value > 0:
            ax.text(
                j,
                i,
                int(value),
                ha="center",
                va="center",
                fontsize=7,
                color="black"
            )

plt.tight_layout()

plt.savefig(
    "H1N1_NA_heatmap_top20_publication.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Final publication heatmap created")
