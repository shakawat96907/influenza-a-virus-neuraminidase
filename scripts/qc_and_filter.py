#!/usr/bin/env python3
"""
Quality control filtering for retrieved CDS nucleotide and protein sequences:
- Protein length within 464-474 amino acids
- Ambiguous bases <= 1%
"""
import pandas as pd
from Bio import SeqIO

def run_qc(in_nuc, in_prot, out_nuc, out_prot, out_map, out_removed):
    nuc_dict = {r.id.split("."): r for r in SeqIO.parse(in_nuc, "fasta")}
    prot_records = list(SeqIO.parse(in_prot, "fasta"))

    kept_nuc, kept_prot, map_rows, removed_rows = [], [], [], []
    for p_rec in prot_records:
        acc = p_rec.id.split(".")
        n_rec = nuc_dict.get(acc)
        if not n_rec:
            continue
        plen = len(p_rec.seq)
        nseq = str(n_rec.seq).upper()
        ambig_pct = (sum(1 for b in nseq if b not in "ACGT") / len(nseq)) * 100

        if not (464 <= plen <= 474):
            removed_rows.append({"Accession": acc, "Reason": "protein_length_out_of_range_464_474"})
            continue
        if ambig_pct > 1.0:
            removed_rows.append({"Accession": acc, "Reason": "ambiguous_bases_gt_1pct"})
            continue

        p_rec.id = acc
        p_rec.description = ""
        n_rec.id = acc
        n_rec.description = ""
        kept_nuc.append(n_rec)
        kept_prot.append(p_rec)
        map_rows.append({"Header": acc, "Accession": acc, "Prot_Len": plen, "Nuc_Len": len(nseq)})

    SeqIO.write(kept_nuc, out_nuc, "fasta")
    SeqIO.write(kept_prot, out_prot, "fasta")
    pd.DataFrame(map_rows).to_csv(out_map, index=False)
    pd.DataFrame(removed_rows).to_csv(out_removed, index=False)
    print(f"total records: {len(prot_records)} | kept: {len(kept_prot)} | removed: {len(removed_rows)}")

if __name__ == "__main__":
    run_qc(
        "raw/new/all_cds_nuc.fasta",
        "raw/new/all_cds_prot.fasta",
        "data/clean_nuc.fasta",
        "data/clean_prot.fasta",
        "data/header_to_accession.csv",
        "results/qc_removed.csv"
    )
