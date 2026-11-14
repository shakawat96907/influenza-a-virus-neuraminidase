from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt

alignment_file = "H1N1_NA_Asia_aligned.fasta"

resistance_markers = {
    275: ("H", ["Y"], "H275Y", "Oseltamivir resistance"),
    223: ("I", ["R", "V", "K"], "I223R/V/K", "Reduced NAI susceptibility"),
    247: ("S", ["N"], "S247N", "Reduced NAI susceptibility"),
    119: ("E", ["V"], "E119V", "Reduced NAI susceptibility"),
    136: ("Q", ["K"], "Q136K", "Reduced NAI susceptibility"),
    292: ("R", ["K"], "R292K", "Reduced NAI susceptibility"),
    295: ("N", ["S"], "N295S", "Reduced NAI susceptibility")
}

print("Loading alignment...")
alignment = AlignIO.read(alignment_file, "fasta")
total_sequences = len(alignment)
print("Total sequences:", total_sequences)

rows = []

for pos, (ref_aa, resistant_aas, marker, relevance) in resistance_markers.items():
    count = 0
    observed = []

    for rec in alignment:
        seq = str(rec.seq)
        if pos - 1 >= len(seq):
            continue

        aa = seq[pos - 1]

        if aa in resistant_aas:
            count += 1
            observed.append(aa)

    freq = (count / total_sequences) * 100

    rows.append({
        "Marker": marker,
        "Position": pos,
        "Expected_AA": ref_aa,
        "Resistance_AA": "/".join(resistant_aas),
        "Detected_Count": count,
        "Total_Sequences": total_sequences,
        "Frequency_percent": freq,
        "Drug_Relevance": relevance
    })

df = pd.DataFrame(rows)
df.to_csv("drug_resistance_barplot_summary.csv", index=False)

print("\nDrug resistance summary:")
print(df)

plt.figure(figsize=(11, 7))

bars = plt.bar(
    df["Marker"],
    df["Frequency_percent"]
)

plt.title(
    "Prevalence of Known Neuraminidase Inhibitor Resistance Markers\nAsian H1N1 NA Sequences",
    fontsize=15,
    fontweight="bold"
)

plt.xlabel("Resistance-associated mutation", fontsize=12, fontweight="bold")
plt.ylabel("Frequency among Asian H1N1 NA sequences (%)", fontsize=12, fontweight="bold")
plt.xticks(rotation=35, ha="right")

for bar, value, count in zip(bars, df["Frequency_percent"], df["Detected_Count"]):
    plt.text(
        bar.get_x() + bar.get_width()/2,
        value + 0.1,
        f"{value:.2f}%\n(n={count})",
        ha="center",
        va="bottom",
        fontsize=9
    )

plt.grid(axis="y", linestyle="--", alpha=0.35)
plt.tight_layout()

plt.savefig("drug_resistance_barplot_H1N1_NA_Asia.png", dpi=600, bbox_inches="tight")

print("\nSaved:")
print("drug_resistance_barplot_summary.csv")
print("drug_resistance_barplot_H1N1_NA_Asia.png")

plt.show()
