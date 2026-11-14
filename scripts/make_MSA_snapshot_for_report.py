from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# =========================
# Input / output
# =========================
msa_file = "H1N1_2015-2025_aligned.fa"
out_dir = Path("msa_snapshot")
out_dir.mkdir(exist_ok=True)

# Important residues from your project
important_positions = [
    15, 52, 75, 188,
    247, 275, 314, 346,
    365, 369, 386, 389,
    416, 432, 450, 453, 454
]

window = 2          # each residue ±2 positions
n_sequences = 20    # representative sequences for snapshot

# =========================
# Read FASTA
# =========================
def read_fasta(path):
    records = []
    header = None
    seq = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)

    if header is not None:
        records.append((header, "".join(seq)))

    return records

records = read_fasta(msa_file)

if len(records) == 0:
    raise ValueError("No sequences found.")

print(f"[INFO] Total sequences loaded: {len(records)}")
print(f"[INFO] Alignment length: {len(records[0][1])}")

# =========================
# Select representative sequences
# evenly spaced from whole dataset
# =========================
indices = np.linspace(0, len(records) - 1, n_sequences, dtype=int)
selected_records = [records[i] for i in indices]

# =========================
# Build selected columns
# =========================
cols = []
for p in important_positions:
    for c in range(p - window, p + window + 1):
        if 1 <= c <= len(records[0][1]):
            cols.append(c)

cols = sorted(set(cols))

# =========================
# Amino acid color map
# =========================
aa_groups = {
    "hydrophobic": list("AVLIMFWY"),
    "polar": list("STNQ"),
    "positive": list("KRH"),
    "negative": list("DE"),
    "special": list("CGP"),
    "gap": ["-", "."]
}

colors = {
    "hydrophobic": "#9ecae1",
    "polar": "#c7e9c0",
    "positive": "#fdd0a2",
    "negative": "#fdae6b",
    "special": "#dadaeb",
    "gap": "#f0f0f0",
    "other": "#ffffff"
}

def aa_color(aa):
    aa = aa.upper()
    for group, aas in aa_groups.items():
        if aa in aas:
            return colors[group]
    return colors["other"]

# =========================
# Plot MSA snapshot
# =========================
fig_height = max(6, n_sequences * 0.35)
fig_width = 16

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

for i, (header, seq) in enumerate(selected_records):
    short_name = header[:18]

    for j, pos in enumerate(cols):
        aa = seq[pos - 1]
        rect_color = aa_color(aa)

        ax.add_patch(
            plt.Rectangle(
                (j, i),
                1,
                1,
                facecolor=rect_color,
                edgecolor="white",
                linewidth=0.4
            )
        )

        ax.text(
            j + 0.5,
            i + 0.5,
            aa,
            ha="center",
            va="center",
            fontsize=7,
            color="black"
        )

# Highlight important residue columns
for p in important_positions:
    if p in cols:
        j = cols.index(p)
        ax.add_patch(
            plt.Rectangle(
                (j, -0.05),
                1,
                n_sequences + 0.1,
                fill=False,
                edgecolor="red",
                linewidth=1.2
            )
        )

# Axis settings
ax.set_xlim(0, len(cols))
ax.set_ylim(0, n_sequences)
ax.invert_yaxis()

ax.set_xticks(np.arange(len(cols)) + 0.5)
ax.set_xticklabels(cols, rotation=90, fontsize=7)

ax.set_yticks(np.arange(n_sequences) + 0.5)
ax.set_yticklabels([r[0][:18] for r in selected_records], fontsize=7)

ax.set_xlabel("Selected NA alignment positions", fontsize=11)
ax.set_ylabel("Representative Asian H1N1 NA sequences", fontsize=11)

ax.set_title(
    "Representative MSA snapshot of selected H1N1 neuraminidase surveillance residues",
    fontsize=13,
    fontweight="bold",
    pad=12
)

ax.tick_params(length=0)

# Remove border
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

png_out = out_dir / "Figure_MSA_snapshot_selected_NA_residues.png"
pdf_out = out_dir / "Figure_MSA_snapshot_selected_NA_residues.pdf"

plt.savefig(png_out, dpi=600, bbox_inches="tight", facecolor="white")
plt.savefig(pdf_out, bbox_inches="tight", facecolor="white")

print("[DONE] Saved:")
print(png_out)
print(pdf_out)
