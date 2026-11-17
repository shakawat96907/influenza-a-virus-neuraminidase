#!/usr/bin/env python3
import numpy as np
import pandas as pd
from PIL import Image, ImageChops
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def main():
    # 1. Load exact verified Supplementary Table 15
    df = pd.read_csv("results/tables/Supplementary Table 15 - (season_matched_vaccine_divergence_469).csv")
    n_leg_gt50 = int((df["Legacy_Ref_Div_Global_Pct"] >= 50.0).sum())
    n_sm_gt50 = int((df["SeasonMatched_Div_Global_Pct"] >= 50.0).sum())

    # Render Clean Figure 9 (1000 DPI TIFF + 1000 DPI PNG)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14.5, 9.0), gridspec_kw={"height_ratios": [1.05, 1.0], "hspace": 0.48})

    l1, = ax1.plot(df["N1_Position"], df["Legacy_Ref_Div_Global_Pct"], color="#95A5A6", lw=1.0, alpha=0.80,
                   label=f"Misassigned Legacy Reference QFQ15359.1 ({n_leg_gt50} sites ≥50%)")
    l2, = ax1.plot(df["N1_Position"], df["SeasonMatched_Div_Global_Pct"], color="#C0392B", lw=1.45, alpha=0.95,
                   label=f"Season-Matched WHO H1N1pdm09 Vaccines ({n_sm_gt50} sites ≥50%)")

    ax1.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Divergence from Vaccine Strain (%)", fontsize=10, fontweight="bold")
    ax1.set_title("A   Elimination of Non-pdm09 Reference Artifacts via Season-Matched WHO Vaccine Comparison (N = 1,449)",
                  loc="left", fontsize=10.6, fontweight="bold", pad=30)
    ax1.set_xlim(0, 474)
    ax1.set_ylim(-3, 108)
    ax1.legend(handles=[l2, l1], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=8.4, frameon=False)
    ax1.grid(True, linestyle=":", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    top_b = df.sort_values("SeasonMatched_Div_Global_Pct", ascending=False).head(15)
    x = np.arange(len(top_b))
    w = 0.36
    b_a = ax2.bar(x - w/2, top_b["SeasonMatched_Div_Asia_Pct"], width=w, color="#C0392B", label="Asian Cohort (n = 750)")
    b_na = ax2.bar(x + w/2, top_b["SeasonMatched_Div_NonAsia_Pct"], width=w, color="#1F618D", label="Non-Asian Control (n = 699)")

    for bars in [b_a, b_na]:
        for bar in bars:
            h = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2, h + 0.6, f"{h:.1f}", ha="center", va="bottom", fontsize=7.0, fontweight="bold")

    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{r['Consensus_AA']}{int(r['N1_Position'])}" for _, r in top_b.iterrows()], fontsize=9, fontweight="bold")
    ax2.set_xlabel("Top Season-Matched Vaccine-Divergent N1 Positions", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Season-Matched Divergence (%)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Regional Comparison of Top Season-Matched WHO Vaccine-Divergent Sites (2015–2025)",
                  loc="left", fontsize=10.6, fontweight="bold", pad=30)
    ax2.set_ylim(0, 46)
    ax2.legend(handles=[b_a, b_na], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=8.4, frameon=False)
    ax2.grid(axis="y", linestyle=":", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.97, top=0.90, bottom=0.08)
    fig.savefig("results/figures/Figure 9/Figure9_SeasonMatched_Vaccine_Divergence_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure 9/Figure9_SeasonMatched_Vaccine_Divergence_1000dpi.png", dpi=1000, format="png")
    plt.close(fig)
    print("1/2: Saved clean Figure 9 (1000 DPI TIFF & PNG) in results/figures/Figure 9/")

    # 2. Render Clean Figure 10 (1000 DPI TIFF + 1000 DPI PNG) from PyMOL render
    im = Image.open("results/structure/pymol_fig10_raw.png").convert("RGB")
    bg = Image.new(im.mode, im.size, (255, 255, 255))
    bbox = ImageChops.difference(im, bg).getbbox()
    if bbox:
        it = iter(bbox)
        x0 = max(0, next(it) - 25)
        y0 = max(0, next(it) - 25)
        x1 = min(im.width, next(it) + 25)
        y1 = min(im.height, next(it) + 25)
        im = im.crop((x0, y0, x1, y1))

    fig2, ax = plt.subplots(figsize=(10.5, 8.2))
    ax.imshow(im)
    ax.axis("off")
    ax.set_title("Structural Localization of Season-Matched WHO Vaccine-Divergent Residues on the N1 Homotetramer (PDB 4B7R)",
                 loc="left", fontsize=10.6, fontweight="bold", pad=52)

    empty_x, empty_y = list(), list()
    leg = [
        Patch(facecolor="#1C6BC7", label="Subunit A (Focus Monomer)"),
        Patch(facecolor="#17A589", label="Subunits B–D (Teal/Coral/Gold)"),
        Line2D(empty_x, empty_y, marker="o", color="w", markerfacecolor="#FF8C00", markersize=8, label="Season-Matched Vaccine-Divergent Sites (≥15%, CPK)"),
        Line2D(empty_x, empty_y, marker="s", color="w", markerfacecolor="#C0392B", markersize=8, label="Conserved Catalytic Pocket (CPK)")
    ]
    ax.legend(handles=leg, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8.4, frameon=False)
    plt.subplots_adjust(left=0.03, right=0.97, top=0.82, bottom=0.04)

    fig2.savefig("results/figures/Figure 10/Figure10_Vaccine_Divergent_3D_Structure_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig2.savefig("results/figures/Figure 10/Figure10_Vaccine_Divergent_3D_Structure_1000dpi.png", dpi=1000, format="png")
    plt.close(fig2)
    print("2/2: Saved clean Figure 10 (1000 DPI TIFF & PNG) in results/figures/Figure 10/")

if __name__ == "__main__":
    main()
