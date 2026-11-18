#!/usr/bin/env python3
"""
Sample year-matched, country-capped global control sequences (cap10, cap20, all)
from deduplicated non-Asian metadata to match annual Asian sample counts.
"""
import pandas as pd

def sample_controls(nonasia_csv, asia_cap20_csv, cap=20, seed=42):
    nonasia = pd.read_csv(nonasia_csv)
    if "Host" in nonasia.columns:
        nonasia = nonasia[nonasia["Host"].str.contains("Homo sapiens", case=False, na=False)]
    asia = pd.read_csv(asia_cap20_csv)
    targets = asia.groupby("Year").size().to_dict()

    capped_pool = nonasia.groupby(["Country", "Year"], group_keys=False).apply(
        lambda g: g.sample(n=min(len(g), cap), random_state=seed)
    )

    selected = []
    for yr, target_n in sorted(targets.items()):
        yr_pool = capped_pool[capped_pool["Year"] == yr]
        n_take = min(len(yr_pool), target_n)
        selected.append(yr_pool.sample(n=n_take, random_state=seed))

    res = pd.concat(selected, ignore_index=True)
    res.to_csv(f"results/control_cap{cap}.csv", index=False)
    res["Accession"].to_csv(f"data/control_human_cap{cap}_accessions.txt", index=False, header=False)
    print(f"cap{cap}: total controls = {len(res)}")

if __name__ == "__main__":
    if pd.io.common.file_exists("results/nonasia_all.csv") and pd.io.common.file_exists("results/asia_selected_cap20.csv"):
        sample_controls("results/nonasia_all.csv", "results/asia_selected_cap20.csv", cap=20)
