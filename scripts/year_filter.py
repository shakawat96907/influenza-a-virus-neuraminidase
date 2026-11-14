from Bio import SeqIO

input_fasta = "H1N1_NA_Asia.fasta"
output_fasta = "H1N1_NA_Asia_2015_2025.fasta"

records = []
for rec in SeqIO.parse(input_fasta, "fasta"):
    header = rec.description
    for year in range(2015, 2026):
        if str(year) in header:
            records.append(rec)
            break

SeqIO.write(records, output_fasta, "fasta")
print("2015–2025 Asia sequences saved:", len(records))
