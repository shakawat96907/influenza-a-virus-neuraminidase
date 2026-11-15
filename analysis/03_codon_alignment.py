#!/usr/bin/env python3
"""
Step 03a (Pipeline Step 7 / T3): Codon-Aware Nucleotide Alignment
Aligns protein sequences via MAFFT (if not already present), trims nucleotide CDS
sequences to exact protein reading frames (removing UTRs/terminal stop codons),
and generates a codon-aware nucleotide alignment using pal2nal.pl.
"""

import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq


def trim_nuc_to_prot(nuc_seq_str, prot_seq_str):
    nuc_clean = nuc_seq_str.upper().replace("-", "").replace(" ", "")
    prot_clean = prot_seq_str.upper().replace("-", "").replace("*", "")
    target_len = len(prot_clean) * 3

    # Check if sequence already matches exact length (with or without terminal stop codon)
    if len(nuc_clean) == target_len:
        return nuc_clean
    if len(nuc_clean) == target_len + 3 and str(Seq(nuc_clean[-3:]).translate()) == "*":
        return nuc_clean[:target_len]

    # Search for exact matching reading frame window if UTRs exist
    best_offset = 0
    best_matches = -1
    max_offset = max(0, len(nuc_clean) - target_len)
    for offset in range(max_offset + 1):
        sub = nuc_clean[offset:offset + target_len]
        if len(sub) < target_len:
            break
        trans = str(Seq(sub).translate())
        matches = sum(a == b for a, b in zip(trans, prot_clean))
        if matches > best_matches:
            best_matches = matches
            best_offset = offset
        if matches == len(prot_clean):
            break

    return nuc_clean[best_offset:best_offset + target_len]


def main():
    os.makedirs("results", exist_ok=True)
    prot_fasta = "data/processed/final_cap20_prot.fasta"
    nuc_fasta = "data/processed/final_cap20_nuc.fasta"
    prot_aln = "results/aln_prot_cap20.fasta"
    trimmed_nuc = "results/trimmed_cds_nuc_cap20.fasta"
    codon_aln = "results/codon_aln_nuc_cap20.fasta"

    # 1. Ensure protein alignment exists
    if not os.path.exists(prot_aln):
        print("[1/3] Running MAFFT protein alignment...")
        with open(prot_aln, "w") as out_f:
            subprocess.run(["mafft", "--auto", "--thread", "4", prot_fasta], stdout=out_f, check=True)
    else:
        print("[1/3] Using existing protein alignment:", prot_aln)

    # 2. Trim nucleotide sequences to exact CDS frame matching protein records
    print("[2/3] Trimming nucleotide sequences to exact protein CDS frames...")
    prot_dict = {r.id: str(r.seq) for r in SeqIO.parse(prot_fasta, "fasta")}
    trimmed_records = []
    for r in SeqIO.parse(nuc_fasta, "fasta"):
        if r.id in prot_dict:
            exact_cds = trim_nuc_to_prot(str(r.seq), prot_dict[r.id])
            r.seq = Seq(exact_cds)
            r.description = ""
            trimmed_records.append(r)

    SeqIO.write(trimmed_records, trimmed_nuc, "fasta")
    print(f"Wrote {len(trimmed_records)} frame-matched CDS records to {trimmed_nuc}")

    # 3. Generate codon alignment using pal2nal.pl (with clean gap-preserving mapping)
    print("[3/3] Generating codon-aware nucleotide alignment via pal2nal.pl...")
    res = subprocess.run(
        ["pal2nal.pl", prot_aln, trimmed_nuc, "-output", "fasta"],
        capture_output=True, text=True
    )

    if res.returncode == 0 and len(res.stdout.strip()) > 0:
        with open(codon_aln, "w") as out_f:
            out_f.write(res.stdout)
        print("pal2nal.pl completed successfully.")
    else:
        print("pal2nal.pl reported ambiguous/stop codons; applying direct frame-guided codon projection...")
        nuc_map = {r.id: str(r.seq) for r in trimmed_records}
        with open(codon_aln, "w") as out_f:
            for rec in SeqIO.parse(prot_aln, "fasta"):
                nseq = nuc_map[rec.id]
                c_idx = 0
                codons = []
                for aa in str(rec.seq):
                    if aa == "-":
                        codons.append("---")
                    else:
                        codons.append(nseq[c_idx:c_idx + 3])
                        c_idx += 3
                out_f.write(f">{rec.id}\n{''.join(codons)}\n")

    subprocess.run(["seqkit", "stats", codon_aln], check=True)


if __name__ == "__main__":
    main()
