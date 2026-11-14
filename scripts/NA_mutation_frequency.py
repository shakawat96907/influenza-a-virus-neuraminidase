from Bio import AlignIO
import pandas as pd

alignment = AlignIO.read("H1N1-2015-2025_aligned_trimmed.fa", "fasta")

reference = alignment[0].seq
total_seq = len(alignment)

mutation_data = []

for pos in range(len(reference)):
    ref_aa = reference[pos]
    if ref_aa == "-":
        continue

    column = [rec.seq[pos] for rec in alignment if rec.seq[pos] != "-"]
    aa_count = {}

    for aa in column:
        aa_count[aa] = aa_count.get(aa, 0) + 1

    for aa, count in aa_count.items():
        if aa != ref_aa:
            mutation_data.append({
                "Position": pos + 1,
                "Reference_AA": ref_aa,
                "Mutant_AA": aa,
                "Count": count,
                "Frequency_%": round((count / total_seq) * 100, 2)
            })

df = pd.DataFrame(mutation_data)
df.to_csv("H1N1_NA_mutation_frequency.csv", index=False)

print("Mutation analysis completed!")
print(df.head())
