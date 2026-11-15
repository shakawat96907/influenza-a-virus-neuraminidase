#!/usr/bin/env python3
import os
import subprocess
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO

def main():
    os.makedirs("results/treetime_clock", exist_ok=True)
    os.makedirs("results/gard", exist_ok=True)

    meta = pd.read_csv("data/processed/final_cap20_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col = col_map.get("accession", list(meta.columns).pop(0))
    date_col = col_map.get("collection_date") or col_map.get("date")
    year_col = col_map.get("year")

    dates_out = []
    for _, row in meta.iterrows():
        acc = str(row[acc_col]).strip()
        raw_d = str(row[date_col]).strip() if date_col and pd.notna(row[date_col]) else ""
        if len(raw_d) >= 10 and "-" in raw_d:
            try:
                yr, mo, dy = map(int, raw_d[:10].split("-"))
                dec = yr + (mo - 1) / 12.0 + (dy - 1) / 365.25
                dates_out.append({"name": acc, "date": f"{dec:.4f}"})
                continue
            except ValueError:
                pass
        if len(raw_d) >= 7 and "-" in raw_d:
            try:
                yr, mo = map(int, raw_d[:7].split("-"))
                dec = yr + (mo - 0.5) / 12.0
                dates_out.append({"name": acc, "date": f"{dec:.4f}"})
                continue
            except ValueError:
                pass
        yr_val = str(row[year_col]).strip() if year_col and pd.notna(row[year_col]) else raw_d[:4]
        dates_out.append({"name": acc, "date": f"{float(yr_val) + 0.5:.4f}"})

    dates_csv = "results/treetime_clock/dates.csv"
    pd.DataFrame(dates_out).to_csv(dates_csv, index=False)
    print(f"[1/3] Exported TreeTime dates for {len(dates_out)} sequences to {dates_csv}")

    print("[2/3] Running TreeTime clock regression & outlier check...")
    cmd = [
        "treetime", "clock",
        "--tree", "results/trees/tree_nuc_cap20_rooted.treefile",
        "--aln", "results/codon_aln_nuc_cap20.fasta",
        "--dates", dates_csv,
        "--clock-filter", "3.0",
        "--outdir", "results/treetime_clock"
    ]
    subprocess.run(cmd, check=True)

    tree = Phylo.read("results/trees/tree_nuc_cap20_rooted.treefile", "newick")
    tips = tree.get_terminals()
    blens = np.array([t.branch_length or 0.0 for t in tips])
    max_tip = max(tips, key=lambda t: t.branch_length or 0.0)
    print("\n=== Terminal Branch Length Audit (1,456 tips) ===")
    print(f"Median terminal branch length : {np.median(blens):.6f}")
    print(f"99th percentile branch length : {np.percentile(blens, 99):.6f}")
    print(f"Max terminal branch length    : {max_tip.branch_length:.6f} (Tip: {max_tip.name})")

    print("\n[3/3] Preparing unique sequences for HyPhy GARD recombination test...")
    seen_seqs = {}
    unique_recs = []
    for rec in SeqIO.parse("results/codon_aln_nuc_cap20.fasta", "fasta"):
        s = str(rec.seq)
        if s not in seen_seqs:
            seen_seqs[s] = rec.id
            unique_recs.append(rec)
    uniq_fasta = "results/gard/unique_codon_aln.fasta"
    SeqIO.write(unique_recs, uniq_fasta, "fasta")
    print(f"Unique nucleotide alignments for GARD: {len(unique_recs)} (from 1,456 total)")

if __name__ == "__main__":
    main()
