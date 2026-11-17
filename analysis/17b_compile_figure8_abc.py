#!/usr/bin/env python3
import pandas as pd
from PIL import Image, ImageChops
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def trim_whitespace(im, pad=62):
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        vals = list(bbox)
        x0 = max(0, vals.pop(0) - pad)
        y0 = max(0, vals.pop(0) - pad)
        x1 = min(im.width, vals.pop(0) + pad)
        y1 = min(im.height, vals.pop(0) + pad)
        return im.crop((x0, y0, x1, y1))
    return im

def main():
    df = pd.read_csv("results/tables/Supplementary Table 14 - (tetramer_vs_monomer_SASA_RSA_4B7R).csv")
    im_3d = trim_whitespace(Image.open("results/structure/pymol_8C_raw.png").convert("RGB"), pad=62)

    # Single horizontal row: 1 row x 3 columns (Panels A, B, C side-by-side)
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(18.5, 5.8),
        gridspec_kw={"width_ratios": [1.0, 1.15, 1.15], "wspace": 0.25}
    )

    # Panel A: Monomer vs Homotetramer SASA
    non_int = df[df["Delta_SASA_Buried"] < 10.0]
    yes_int = df[df["Delta_SASA_Buried"] >= 10.0]
    s1 = ax1.scatter(non_int["Monomer_SASA_A2"], non_int["Tetramer_SASA_A2"], s=16, color="#2E86C1", alpha=0.72, label=f"Non-interfacial (n={len(non_int)})")
    s2 = ax1.scatter(yes_int["Monomer_SASA_A2"], yes_int["Tetramer_SASA_A2"], s=25, color="#C0392B", alpha=0.88, edgecolors="#111111", lw=0.4, label=f"Interfacial (ΔSASA≥10 Å², n={len(yes_int)})")
    diag_pts = list(range(0, 225, 25))
    ax1.plot(diag_pts, diag_pts, color="#333333", linestyle="--", lw=1.0)

    for _, r in yes_int.nlargest(4, "Delta_SASA_Buried").iterrows():
        ax1.annotate(
            r["Residue_4B7R"],
            xy=(r["Monomer_SASA_A2"], r["Tetramer_SASA_A2"]),
            xytext=(r["Monomer_SASA_A2"] + 10, max(3, r["Tetramer_SASA_A2"] - 5)),
            fontsize=7.8, fontweight="bold", color="#78281F",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#555555")
        )

    ax1.set_xlabel("Isolated Monomer SASA (Å², Chain A)", fontsize=9.5, fontweight="bold")
    ax1.set_ylabel("Biological Tetramer SASA (Å², 4B7R)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A   Monomer vs. Tetramer SASA (4B7R)", loc="left", fontsize=10.5, fontweight="bold", pad=62)
    ax1.set_xlim(-5, 225)
    ax1.set_ylim(-5, 225)
    ax1.legend(handles=[s1, s2], loc="lower left", bbox_to_anchor=(-0.02, 1.01), ncol=2, fontsize=7.5, frameon=False)
    ax1.grid(True, linestyle=":", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Biological Tetramer RSA Profile
    exp_hi = df[df["Tetramer_RSA"] >= 0.25]
    exp_md = df[(df["Tetramer_RSA"] >= 0.05) & (df["Tetramer_RSA"] < 0.25)]
    exp_bu = df[df["Tetramer_RSA"] < 0.05]

    b1 = ax2.scatter(exp_bu["N1_Position"], exp_bu["Tetramer_RSA"], s=13, color="#95A5A6", alpha=0.65, label=f"Buried (<0.05, n={len(exp_bu)})")
    b2 = ax2.scatter(exp_md["N1_Position"], exp_md["Tetramer_RSA"], s=18, color="#F39C12", alpha=0.82, label=f"Moderate (0.05–0.25, n={len(exp_md)})")
    b3 = ax2.scatter(exp_hi["N1_Position"], exp_hi["Tetramer_RSA"], s=25, color="#27AE60", alpha=0.88, label=f"Exposed (≥0.25, n={len(exp_hi)})")
    ax2.axhline(0.25, color="#1E8449", linestyle="--", lw=0.9)
    ax2.axhline(0.05, color="#7F8C8D", linestyle=":", lw=0.9)

    label_pos = {200, 222, 275, 339, 389, 416, 453, 469}
    for idx_k, (_, r) in enumerate(df[df["N1_Position"].isin(label_pos)].iterrows()):
        dy = 0.08 if idx_k % 2 == 0 else -0.10
        ax2.annotate(
            r["Residue_4B7R"],
            xy=(r["N1_Position"], r["Tetramer_RSA"]),
            xytext=(r["N1_Position"], min(0.95, max(0.03, r["Tetramer_RSA"] + dy))),
            fontsize=7.8, fontweight="bold", ha="center", color="#1A2530",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#555555")
        )

    ax2.set_xlabel("Canonical N1 Head-Domain Position (83–469)", fontsize=9.5, fontweight="bold")
    ax2.set_ylabel("Tetramer Relative Solvent Accessibility (RSA)", fontsize=9.5, fontweight="bold")
    ax2.set_title("B   Biological Tetramer RSA Exposure Profile", loc="left", fontsize=10.5, fontweight="bold", pad=62)
    ax2.set_xlim(78, 475)
    ax2.set_ylim(-0.04, 1.05)
    ax2.legend(handles=[b3, b2, b1], loc="lower left", bbox_to_anchor=(-0.02, 1.01), ncol=3, fontsize=7.3, frameon=False)
    ax2.grid(True, linestyle=":", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Panel C: PyMOL 3D Biological Homotetramer
    ax3.imshow(im_3d)
    ax3.axis("off")
    ax3.set_title("C   3D Homotetramer Architecture (PDB 4B7R)", loc="left", fontsize=10.5, fontweight="bold", pad=62)

    empty_x = list()
    empty_y = list()
    leg_c = [
        Patch(facecolor="#1C6BC7", edgecolor="none", label="Subunit A"),
        Patch(facecolor="#17A589", edgecolor="none", label="Subunits B–D (Teal/Coral/Gold)"),
        Line2D(empty_x, empty_y, marker="o", color="w", markerfacecolor="#FF8C00", markersize=7, label="Head Variants (CPK)"),
        Line2D(empty_x, empty_y, marker="s", color="w", markerfacecolor="#C0392B", markersize=7, label="Catalytic Site (CPK)")
    ]
    ax3.legend(handles=leg_c, loc="lower left", bbox_to_anchor=(-0.02, 1.01), ncol=2, fontsize=7.4, frameon=False)

    plt.subplots_adjust(left=0.045, right=0.985, top=0.74, bottom=0.13)
    fig.savefig("results/figures/Figure8_Composite_Tetramer_Structure_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure8_Composite_Tetramer_Structure_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("Saved 1-row 3-panel Figure 8 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
