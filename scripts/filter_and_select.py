#!/usr/bin/env python3
"""
Filter Asian H1N1 NA metadata (2015-2025, complete CDS 1380-1500 nt, human host),
deduplicate by Isolate name, and perform country-year capped sampling (cap10, cap20, all).
"""
import pandas as pd

def filter_asian_metadata(meta_csv, out_prefix="data/asia_selected"):
    df = pd.read_csv(meta_csv)
    df["Year"] = pd.to_datetime(df["Collection_Date"], errors="coerce").dt.year
    df = df[df["Year"].between(2015, 2025)]
    df = df[df["Length"].between(1380, 1500)]
    if "Host" in df.columns:
        df = df[df["Host"].str.contains("Homo sapiens", case=False, na=False)]

    # Deduplicate by Isolate name (keeping one record per isolate)
    has_iso = df["Isolate"].notna() & (df["Isolate"].astype(str).str.strip() != "")
    df_iso = df[has_iso].drop_duplicates(subset=["Isolate"], keep="first")
    df_no_iso = df[~has_iso].drop_duplicates(subset=["Accession"], keep="first")
    dedup = pd.concat([df_iso, df_no_iso], ignore_index=True)

    dedup["Accession"].to_csv(f"{out_prefix}_all_accessions.txt", index=False, header=False)
    for cap in (10, 20):
        sampled = dedup.groupby(["Country", "Year"], group_keys=False).apply(
            lambda g: g.sample(n=min(len(g), cap), random_state=42)
        )
        sampled["Accession"].to_csv(f"{out_prefix}_cap{cap}_accessions.txt", index=False, header=False)
    print(f"Asian deduplicated records: {len(dedup)}")

if __name__ == "__main__":
    if pd.io.common.file_exists("raw/ncbi_asia_metadata.csv"):
        filter_asian_metadata("raw/ncbi_asia_metadata.csv")
