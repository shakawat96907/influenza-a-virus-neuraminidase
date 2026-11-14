import pandas as pd

df = pd.read_csv("H1N1_NA_mutation_frequency.csv")

dom = df.sort_values("Count", ascending=False) \
        .groupby("Position").first().reset_index()

dom.to_csv("H1N1_NA_Asia_dominant_mutations.csv", index=False)
print("Dominant mutations saved")
