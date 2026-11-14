from Bio import SeqIO
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
import os

# =========================
# INPUT / SETTINGS
# =========================
alignment_file = "H1N1_2015-2025_aligned.fa"   # alignment file name
output_dir = "appendix_msa_figures"

# total representative sequences to show
n_total_sequences = 30

# split into 2 figures
n_per_panel = 15

# important selected residues (1-based numbering)
selected_positions = [13, 15, 51, 52, 74, 75, 77, 81,
                      188, 314, 346, 365, 369, 386, 389,
                      416, 432, 450, 453, 454]

main_title = "Representative MSA snapshot of selected H1N1 neuraminidase residues"

# =========================
# AMINO ACID COLOR MAP
# =========================
aa_colors = {
    'A': '#d9f0d3', 'V': '#ccebc5', 'L': '#b3e2cd', 'I': '#ccebc5',
    'M': '#ccebc5', 'F': '#fddaec', 'W': '#f1eef6', 'Y': '#fff2ae',
    'S': '#c6dbef', 'T': '#bdd7e7', 'N': '#dadaeb', 'Q': '#dadaeb',
    'K': '#fee391', 'R': '#fec44f', 'H': '#fdd0a2',
    'D': '#fcbba1', 'E': '#fcae91',
    'G': '#e5f5e0', 'P': '#f2f0f7', 'C': '#fdd49e',
    '-': '#ffffff', 'X': '#f0f0f0'
}

def get_color(aa):
    return aa_colors.get(aa.upper(), '#f7f7f7')

# =========================
# LOAD ALIGNMENT
# =========================
records = list(SeqIO.parse(alignment_file, "fasta"))
if len(records) == 0:
    raise ValueError("No sequences found in alignment file.")

alignment_length = len(records[0].seq)

# check all aligned lengths are same
for rec in records:
    if len(rec.seq) != alignment_length:
        raise ValueError("Sequences are not equally aligned. Check your alignment file.")

print(f"[INFO] Total aligned sequences found: {len(records)}")
print(f"[INFO] Alignment length: {alignment_length}")

# =========================
# SELECT 30 REPRESENTATIVE SEQUENCES
# evenly spaced across whole dataset
# =========================
if len(records) <= n_total_sequences:
    selected_records = records
else:
    step = (len(records) - 1) / (n_total_sequences - 1)
    indices = sorted(set([round(i * step) for i in range(n_total_sequences)]))
    selected_records = [records[i] for i in indices]

# if rounding produced slightly fewer, fill manually
while len(selected_records) < n_total_sequences and len(selected_records) < len(records):
    for rec in records:
        if rec not in selected_records:
            selected_records.append(rec)
        if len(selected_records) == n_total_sequences:
            break

selected_records = selected_records[:n_total_sequences]

print(f"[INFO] Selected representative sequences: {len(selected_records)}")

# =========================
# SHORTEN NAMES
# =========================
def short_name(header):
    h = header.strip()
    if len(h) > 18:
        return h[:18]
    return h

seq_names = [short_name(r.id) for r in selected_records]

# =========================
# EXTRACT SELECTED POSITIONS
# =========================
# convert 1-based positions to 0-based
pos_idx = [p - 1 for p in selected_positions]

def extract_matrix(record_subset):
    matrix = []
    names = []
    for rec in record_subset:
        row = []
        seq = str(rec.seq)
        for idx in pos_idx:
            if idx < len(seq):
                row.append(seq[idx])
            else:
                row.append('-')
        matrix.append(row)
        names.append(short_name(rec.id))
    return matrix, names

# =========================
# PLOTTING FUNCTION
# =========================
def plot_panel(matrix, names, panel_no, output_png, output_pdf):
    n_rows = len(matrix)
    n_cols = len(matrix[0])

    # A4 landscape-ish size
    fig_w = 11.69
    fig_h = 8.27

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=300)

    # margins
    left_margin = 3.2
    bottom_margin = 1.4
    cell_w = 0.55
    cell_h = 0.42

    # draw cells
    for i in range(n_rows):
        for j in range(n_cols):
            aa = matrix[i][j]
            x = left_margin + j * cell_w
            y = bottom_margin + (n_rows - 1 - i) * cell_h

            rect = Rectangle((x, y), cell_w, cell_h,
                             facecolor=get_color(aa),
                             edgecolor='gray', linewidth=0.35)
            ax.add_patch(rect)

            ax.text(x + cell_w/2, y + cell_h/2, aa,
                    ha='center', va='center',
                    fontsize=6.5, color='black')

    # row labels
    for i, nm in enumerate(names):
        y = bottom_margin + (n_rows - 1 - i) * cell_h + cell_h/2
        ax.text(left_margin - 0.15, y, nm,
                ha='right', va='center', fontsize=7.5)

    # column labels
    for j, p in enumerate(selected_positions):
        x = left_margin + j * cell_w + cell_w/2
        ax.text(x, bottom_margin - 0.08, str(p),
                ha='center', va='top', fontsize=7, rotation=90)

    # axis labels
    ax.text(left_margin + (n_cols * cell_w)/2,
            bottom_margin - 0.55,
            "Selected NA alignment positions",
            ha='center', va='top', fontsize=10, fontweight='bold')

    ax.text(left_margin - 1.35,
            bottom_margin + (n_rows * cell_h)/2,
            "Representative Asian H1N1 NA sequences",
            ha='center', va='center', rotation=90,
            fontsize=10, fontweight='bold')

    # title
    ax.text(left_margin + (n_cols * cell_w)/2,
            bottom_margin + n_rows * cell_h + 0.55,
            f"{main_title} (Panel {panel_no})",
            ha='center', va='bottom',
            fontsize=13, fontweight='bold')

    # limits / clean
    ax.set_xlim(0, left_margin + n_cols * cell_w + 0.6)
    ax.set_ylim(0, bottom_margin + n_rows * cell_h + 1.0)
    ax.axis('off')

    plt.tight_layout()
    plt.savefig(output_png, dpi=600, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()

    print(f"[DONE] Saved: {output_png}")
    print(f"[DONE] Saved: {output_pdf}")

# =========================
# SPLIT INTO 2 PANELS
# =========================
os.makedirs(output_dir, exist_ok=True)

panel1_records = selected_records[:n_per_panel]
panel2_records = selected_records[n_per_panel:n_per_panel*2]

matrix1, names1 = extract_matrix(panel1_records)
matrix2, names2 = extract_matrix(panel2_records)

plot_panel(
    matrix1, names1, 1,
    os.path.join(output_dir, "Appendix_MSA_panel1_15seq.png"),
    os.path.join(output_dir, "Appendix_MSA_panel1_15seq.pdf")
)

plot_panel(
    matrix2, names2, 2,
    os.path.join(output_dir, "Appendix_MSA_panel2_15seq.png"),
    os.path.join(output_dir, "Appendix_MSA_panel2_15seq.pdf")
)

print("\n[ALL DONE] Two appendix-friendly MSA figures created.")
