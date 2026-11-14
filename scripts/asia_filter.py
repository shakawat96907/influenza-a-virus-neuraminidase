import sys
import pandas as pd
from Bio import SeqIO

# Command line arguments
fasta_file = sys.argv[1]
metadata_file = sys.argv[2]
output_fasta = sys.argv[3]

# Read metadata (CSV file)
meta = pd.read_csv(metadata_file)

# Asia countries list (edit/add if needed)
asia_countries = [
    "India", "Japan", "China", "Bangladesh", "Thailand",
    "Vietnam", "South Korea", "Pakistan", "Sri Lanka",
    "Nepal", "Indonesia", "Philippines", "Malaysia"
]

# Filter metadata for Asia
asia_meta = meta[meta["country"].isin(asia_countries)]

# Get accession IDs of Asia sequences
asia_accessions = set(asia_meta["accession"].astype(str))

print("Total metadata:", len(meta))
print("Asia metadata:", len(asia_meta))

# Read FASTA and filter
asia_sequences = []

for record in SeqIO.parse(fasta_file, "fasta"):
    # accession is assumed to be first part of FASTA header
    accession = record.id.split("|")[0]
    if accession in asia_accessions:
        asia_sequences.append(record)

# Write output FASTA
SeqIO.write(asia_sequences, output_fasta, "fasta")

print(f"Asia NA sequences saved: {len(asia_sequences)}")
