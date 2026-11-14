#!/usr/bin/env python3
"""
Combine regional NCBI Virus metadata tables (Africa, Europe, North America, Oceania, South America),
filter by year (2015-2025) and complete CDS length (1380-1500 nt), and deduplicate by Isolate.
"""
import os
import pandas as pd

REGIONS = {
    "africa": "raw/ncbi_africa_metadata.csv",
    "europe": "raw/ncbi_europe_metadata.csv",
    "northamerica": "raw/ncbi_northamerica_metadata.csv",
    "oceania": "raw/ncbi_oceania_metadata.csv",
    "southamerica": "raw/ncbi_southamerica_metadata.csv",
}

def main():
    os.makedirs("results", exist_ok=True)
    frames = []
    for reg, path in REGIONS.items():
        if os.path.exists(path):
            df = pd.read_csv(path)
            df["region"] = reg
            frames.append(df)

    combined = pd.concat(frames, ignore_index=True)
    combined["Year"] = pd.to_datetime(combined["Collection_Date"], errors="coerce").dt.year
    combined = combined[combined["Year"].between(2015, 2025)]
    combined = combined[combined["Length"].between(1380, 1500)]

    has_iso = combined["Isolate"].notna() & (combined["Isolate"].astype(str).str.strip() != "")
    dedup = pd.concat([
        combined[has_iso].drop_duplicates(subset=["Isolate"], keep="first"),
        combined[~has_iso].drop_duplicates(subset=["Accession"], keep="first")
    ], ignore_index=True)

    dedup.to_csv("results/nonasia_all.csv", index=False)
    counts = dedup.groupby(["region", "Year"]).size().reset_index(name="count")
    counts.to_csv("results/nonasia_counts_region_year.csv", index=False)
    print(f"Combined deduplicated non-Asian records: {len(dedup)}")

if __name__ == "__main__":
    main()
