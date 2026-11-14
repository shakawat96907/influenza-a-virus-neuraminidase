import pandas as pd
import matplotlib.pyplot as plt

# Load mutation frequency data
df = pd.read_csv("H1N1_NA_mutation_frequency.csv")

# Total mutation count per position
landscape_df = df.groupby("Position")["Count"].sum().reset_index()

# Sort by position
landscape_df = landscape_df.sort_values("Position")

# Plot mutation landscape
plt.figure(figsize=(12,5))

plt.plot(
    landscape_df["Position"],
    landscape_df["Count"],
    color="#4C72B0",
    linewidth=2
)

plt.xlabel("NA Position", fontsize=12)
plt.ylabel("Mutation Count", fontsize=12)

plt.title(
    "Mutation Landscape of H1N1 Neuraminidase (Asia, 2015–2025)",
    fontsize=14
)

plt.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()

plt.savefig("H1N1_NA_mutation_landscape.png", dpi=300)

plt.show()

print("Mutation landscape plot created")

