import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("H1N1_NA_mutation_frequency.csv")

plt.figure(figsize=(14,6))
plt.scatter(df["Position"], df["Count"], alpha=0.6)
plt.xlabel("NA Position")
plt.ylabel("Mutation Count")
plt.title("H1N1 NA Asia Mutation Frequency")
plt.tight_layout()
plt.savefig("H1N1_NA_Asia_mutation_plot.png", dpi=300)
plt.show()
