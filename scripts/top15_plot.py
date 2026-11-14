import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("H1N1_NA_mutation_frequency.csv")

# Same position-এর সব mutation count যোগ করা
pos_counts = df.groupby("Position")["Count"].sum().reset_index()
pos_counts = pos_counts.sort_values("Position")

# Top 15 highest mutated NA positions
top15 = pos_counts.sort_values("Count", ascending=False).head(15)
top15_positions = set(top15["Position"])

print("Top 15 NA positions:")
print(top15)

top15.to_csv("top15_NA_positions.csv", index=False)

plt.figure(figsize=(28, 9))

bars = plt.bar(
    pos_counts["Position"],
    pos_counts["Count"],
    width=0.8
)

max_y = pos_counts["Count"].max()

# শুধু top 15 bar-এর উপরে NA position number
for bar, pos, count in zip(bars, pos_counts["Position"], pos_counts["Count"]):
    if pos in top15_positions:
        plt.text(
            bar.get_x() + bar.get_width()/2,
            count + max_y * 0.02,
            str(int(pos)),
            ha="center",
            va="bottom",
            fontsize=10,
            rotation=90,
            fontweight="bold"
        )

plt.xlabel("NA Position")
plt.ylabel("Mutation Count")
plt.title("H1N1 NA Mutation Frequencies (All Positions)")

plt.xticks(pos_counts["Position"][::10], rotation=90, fontsize=7)
plt.ylim(0, max_y * 1.18)
plt.grid(axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig("H1N1_NA_all_bar_plot_top15_labeled.png", dpi=600)
plt.savefig("H1N1_NA_all_bar_plot_top15_labeled.pdf")
plt.show()
