#!/usr/bin/env python3
"""
Step 07 (Pipeline Step 8 / T4): Outlier Removal, Recombination Screening (HyPhy GARD),
and Clean Dataset Re-Alignment.
Removes the 7 confirmed non-pdm09 / misannotated / extreme long-branch outliers
(MK502498, OM679430, PV156074, OM736137, PP873819, MN508986, MK226524),
yielding 1,449 clean human A(H1N1)pdm09 NA sequences (750 Asian + 699 Non-Asian).
"""

import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

OUTLIERS = {
    "MK502498": "Misannotated non-N1 subtype / extreme long branch (21.54% identity; branch=8.3057)",
    "OM679430": "Swine-origin variant / pre-2009 spillover (76.12% identity; branch=0.1024)",
    "PV156074": "Pre-2009 / variant spillover (80.81% identity; branch=0.2129)",
    "OM736137": "Divergent non-pdm09 variant (86.14% identity)",
    "PP873819": "Divergent non-pdm09 variant (87.21% identity)",
    "MN508986": "Divergent variant / clock outlier (88.70% identity; branch=0.0625)",
    "MK226524": "High-divergence clock outlier (90.19% identity; residual=+14.16 yr)"
}

def main():
    os.makedirs("data/processed", exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)

    # Save documented outlier audit log (for Table 2 & Editor Response)
    out_df = pd.DataFrame([{"Accession": k, "Exclusion_Reason": v} for k, v in OUTLIERS.items()])
    out_df.to_csv("results/tables/removed_phylogenetic_outliers.csv", index=False)

    # Filter metadata, protein FASTA, and nucleotide FASTA
    meta = pd.read_csv("data/processed/final_cap20_meta.csv")
    acc_col = meta.columns
    clean_meta = meta[~meta[acc_col].astype(str).isin(OUTLIERS)].copy()
    clean_meta.to_csv("data/processed/final_clean_meta.csv", index=False)

    clean_prot = [r for r in SeqIO.parse("data/processed/final_cap20_prot.fasta", "fasta") if r.id not in OUTLIERS]
    clean_nuc = [r for r in SeqIO.parse("results/trimmed_cds_nuc_cap20.fasta", "fasta") if r.id not in OUTLIERS]

    SeqIO.write(clean_prot, "data/processed/final_clean_prot.fasta", "fasta")
    SeqIO.write(clean_nuc, "data/processed/final_clean_nuc.fasta", "fasta")

    grp_col = [c for c in clean_meta.columns if c.lower() == "group"]
    counts = clean_meta[grp_col].value_counts().to_dict()
    print(f"=== Clean Dataset Summary (After T4 Outlier Removal) ===")
    print(f"Total clean sequences: {len(clean_meta)} | Breakdown: {counts}")

    # Re-run MAFFT protein alignment & consensus on clean 1,449 sequences
    print("\n[1/3] Re-running MAFFT on 1,449 clean protein sequences...")
    with open("results/aln_prot_clean.fasta", "w") as out_f:
        subprocess.run(["mafft", "--auto", "--thread", "4", "data/processed/final_clean_prot.fasta"], stdout=out_f, check=True)

    # Re-run pal2nal codon alignment on clean 1,449 sequences
    print("[2/3] Generating clean codon-aware nucleotide alignment via pal2nal.pl...")
    res = subprocess.run(
        ["pal2nal.pl", "results/aln_prot_clean.fasta", "data/processed/final_clean_nuc.fasta", "-output", "fasta"],
        capture_output=True, text=True, check=True
    )
    with open("results/codon_aln_nuc_clean.fasta", "w") as out_f:
        out_f.write(res.stdout)

    subprocess.run(["seqkit", "stats", "results/codon_aln_nuc_clean.fasta"], check=True)

    # Prepare 150-sequence representative subset for fast HyPhy GARD recombination screening
    print("\n[3/3] Running HyPhy GARD recombination screening on representative subset...")
    uniq_recs = []
    seen = set()
    for r in SeqIO.parse("results/codon_aln_nuc_clean.fasta", "fasta"):
        s = str(r.seq)
        if s not in seen:
            seen.add(s)
            uniq_recs.append(r)
    step = max(1, len(uniq_recs) // 120)
    sub_recs = uniq_recs[::step][:120]
    sub_fasta = "results/gard/clean_subset_120.fasta"
    SeqIO.write(sub_recs, sub_fasta, "fasta")

    subprocess.run(
        ["hyphy", "gard", "--alignment", sub_fasta, "--output", "results/gard/gard_clean_results.json"],
        check=True
    )
    print("HyPhy GARD recombination test completed successfully!")

if __name__ == "__main__":
    main()
