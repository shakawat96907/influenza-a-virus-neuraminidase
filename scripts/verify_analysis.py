#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO

def main():
    nuc_ids = [r.id for r in SeqIO.parse("data/processed/final_cap20_nuc.fasta", "fasta")]
    prot_ids = [r.id for r in SeqIO.parse("data/processed/final_cap20_prot.fasta", "fasta")]
    meta = pd.read_csv("data/processed/final_cap20_meta.csv")
    assert set(nuc_ids) == set(prot_ids) == set(meta["Accession"].astype(str))
    assert len(nuc_ids) == len(set(nuc_ids))
    print(f"Verified {len(nuc_ids)} unique H1N1 NA records successfully.")

if __name__ == "__main__":
    main()
