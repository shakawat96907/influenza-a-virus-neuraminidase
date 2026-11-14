#!/usr/bin/env python3
"""
Retrieve complete coding sequences (CDS) for nucleotide and protein records
from NCBI GenBank/Entrez in batches of 100 accessions.
"""
import os
import time
from Bio import Entrez, SeqIO

Entrez.email = "research@sust.edu"

def fetch_cds_batches(acc_file, out_nuc, out_prot, batch_size=100):
    with open(acc_file) as f:
        accessions = [line.strip() for line in f if line.strip()]
    print(f"accessions requested: {len(accessions)}")

    os.makedirs(os.path.dirname(out_nuc), exist_ok=True)
    nuc_records, prot_records = [], []

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        print(f"batch {i // batch_size + 1}/{(len(accessions) + batch_size - 1) // batch_size}: {len(batch)} accessions")
        handle = Entrez.efetch(db="nuccore", id=",".join(batch), rettype="gb", retmode="text")
        for record in SeqIO.parse(handle, "genbank"):
            for feat in record.features:
                if feat.type == "CDS" and "translation" in feat.qualifiers:
                    nuc_seq = feat.extract(record.seq)
                    prot_seq = feat.qualifiers["translation"]
                    rec_nuc = record[:0]
                    rec_nuc.id = record.id.split(".")
                    rec_nuc.description = ""
                    rec_nuc.seq = nuc_seq
                    nuc_records.append(rec_nuc)

                    rec_prot = record[:0]
                    rec_prot.id = record.id.split(".")
                    rec_prot.description = ""
                    rec_prot.seq = prot_seq
                    prot_records.append(rec_prot)
                    break
        handle.close()
        time.sleep(0.4)

    SeqIO.write(nuc_records, out_nuc, "fasta")
    SeqIO.write(prot_records, out_prot, "fasta")
    print(f"wrote {out_nuc}: {len(nuc_records)} nucleotide CDS records")
    print(f"wrote {out_prot}: {len(prot_records)} protein CDS records")

if __name__ == "__main__":
    fetch_cds_batches(
        "data/all_union_accessions.txt",
        "raw/new/all_cds_nuc.fasta",
        "raw/new/all_cds_prot.fasta"
    )
