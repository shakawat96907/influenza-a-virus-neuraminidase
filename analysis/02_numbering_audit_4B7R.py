#!/usr/bin/env python3
import os
import urllib.request
from Bio import SeqIO, Align, PDB
from Bio.SeqUtils import seq1
import pandas as pd

def main():
    os.makedirs("results/tables", exist_ok=True)
    pdb_path = "results/4B7R.pdb"
    if not os.path.exists(pdb_path):
        url = "https://files." + "rcsb.org/download/4B7R.pdb"
        urllib.request.urlretrieve(url, pdb_path)

    cons_rec = next(SeqIO.parse("results/consensus_prot_cap20.fasta", "fasta"))
    cons_aln_seq = str(cons_rec.seq)
    cons_nogap_seq = cons_aln_seq.replace("-", "")

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("4B7R", pdb_path)
    chain_a = structure["A"]

    pdb_residues = [(res.id, seq1(res.resname)) for res in chain_a if res.id == " "]
    pdb_seq = "".join([aa for _, aa in pdb_residues])

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    aln = aligner.align(cons_nogap_seq, pdb_seq)
    aligned_cons, aligned_pdb = aln, aln

    matches = sum(a == b for a, b in zip(aligned_cons, aligned_pdb) if a != "-" and b != "-")
    aligned_len = sum(1 for a, b in zip(aligned_cons, aligned_pdb) if a != "-" and b != "-")
    print(f"Sequence Identity: {(matches / aligned_len) * 100:.2f}% ({matches}/{aligned_len})")

    n1_to_pdb = {}
    n1_idx, pdb_idx = 0, 0
    for a, b in zip(aligned_cons, aligned_pdb):
        if a != "-": n1_idx += 1
        if b != "-": pdb_idx += 1
        if a != "-" and b != "-":
            n1_to_pdb[n1_idx] = pdb_residues[pdb_idx - 1]

    rows = []
    n1_pos = 0
    for aln_pos, char in enumerate(cons_aln_seq, start=1):
        if char != "-":
            n1_pos += 1
            pdb_resnum, pdb_aa = n1_to_pdb.get(n1_pos, (None, "-"))
            n2_approx = (n1_pos - 1) if n1_pos >= 65 else n1_pos
            rows.append({
                "MAFFT_Alignment_Pos": aln_pos,
                "Consensus_AA": char,
                "Canonical_N1_Pos": n1_pos,
                "PDB_4B7R_ResNum": pdb_resnum,
                "PDB_4B7R_AA": pdb_aa,
                "Classical_N2_Pos": n2_approx,
            })
        else:
            rows.append({
                "MAFFT_Alignment_Pos": aln_pos,
                "Consensus_AA": "-",
                "Canonical_N1_Pos": None,
                "PDB_4B7R_ResNum": None,
                "PDB_4B7R_AA": "-",
                "Classical_N2_Pos": None,
            })

    pd.DataFrame(rows).to_csv("results/tables/numbering_audit_table.csv", index=False)

if __name__ == "__main__":
    main()
