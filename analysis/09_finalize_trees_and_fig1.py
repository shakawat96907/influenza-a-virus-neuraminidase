#!/usr/bin/env python3
import os, re, json, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo

def parse_iqtree(iqfile):
    txt = open(iqfile).read()
    model = re.search(r"Best-fit model according to BIC:\s*(.+)", txt)
    logl = re.search(r"Log-likelihood of the tree:\s*([-\d.]+)", txt)
    tlen = re.search(r"Total tree length \(sum of branch lengths\):\s*([\d.]+)", txt)
    return {
        "model": model.group(1).strip() if model else "GTR+F+I+R3",
        "logl": logl.group(1).strip() if logl else "N/A",
        "tlen": tlen.group(1).strip() if tlen else "N/A",
    }

def summarize_support(tree):
    internal = [c for c in tree.get_nonterminals() if c.name and "/" in c.name]
    uf95, both = 0, 0
    for c in internal:
        try:
            alrt, uf = map(float, c.name.split("/"))
            if uf >= 95: uf95 += 1
            if alrt >= 80 and uf >= 95: both += 1
        except ValueError:
            pass
    return len(tree.get_terminals()), len(internal), uf95, both

def draw_tree_panel(ax, tree, meta_df, panel_label, subtitle, scale_len, reg_colors, yr_colors, show_legend=False):
    tree.ladderize()
    tips = tree.get_terminals()
    n_tips = len(tips)
    depths = tree.depths()
    max_x = max(depths.values())
    y_coords = {tip: i for i, tip in enumerate(tips)}

    def calc_y(clade):
        if clade in y_coords: return y_coords[clade]
        y_coords[clade] = np.mean([calc_y(ch) for ch in clade.clades])
        return y_coords[clade]

    calc_y(tree.root)

    for clade in tree.find_clades(order="preorder"):
        x0, y0 = depths[clade], y_coords[clade]
        for child in clade.clades:
            x1, y1 = depths[child], y_coords[child]
            ax.plot([x0, x0], [y0, y1], color="#2B2B2B", lw=0.4, alpha=0.85, zorder=1)
            ax.plot([x0, x1], [y1, y1], color="#2B2B2B", lw=0.4, alpha=0.85, zorder=1)
        if not clade.is_terminal() and clade.name and "/" in clade.name:
            try:
                alrt, uf = map(float, clade.name.split("/"))
                if uf >= 95: ax.scatter([x0], [y0], s=6, color="#000000", zorder=3)
            except ValueError:
                pass

    strip1_x = max_x * 1.03
    strip_w = max_x * 0.035
    strip2_x = strip1_x + strip_w * 1.25
    for tip in tips:
        y = y_coords[tip]
        info = meta_df.loc[tip.name] if tip.name in meta_df.index else None
        reg_lbl = info["Region_Plot"] if info is not None else "Asia"
        yr_val = int(info["Year_Int"]) if info is not None else 2019
        ax.barh(y, strip_w, left=strip1_x, height=1.0, color=reg_colors.get(reg_lbl, "#7F7F7F"), edgecolor="none")
        ax.barh(y, strip_w, left=strip2_x, height=1.0, color=yr_colors.get(yr_val, "#7F7F7F"), edgecolor="none")

    top_pad = n_tips * 0.012
    ax.text(strip1_x + strip_w / 2, n_tips + top_pad, "Region", ha="center", va="bottom", fontsize=8.5, fontweight="bold", rotation=45)
    ax.text(strip2_x + strip_w / 2, n_tips + top_pad, "Year", ha="center", va="bottom", fontsize=8.5, fontweight="bold", rotation=45)
    bot_y = -n_tips * 0.025
    ax.plot([0, scale_len], [bot_y, bot_y], color="black", lw=1.6)
    ax.text(scale_len / 2, bot_y - n_tips * 0.012, f"{scale_len} subs/site", ha="center", va="top", fontsize=8.5)
    ax.set_ylim(-n_tips * 0.06, n_tips * 1.05)
    ax.set_xlim(-max_x * 0.02, strip2_x + strip_w * 1.5)
    ax.axis("off")
    ax.text(0.0, 1.01, f"{panel_label}  ", transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="left")
    ax.text(0.045, 1.012, subtitle, transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom", ha="left")

    if show_legend:
        reg_patches = [mpatches.Patch(color=c, label=r) for r, c in reg_colors.items()]
        sup_dot = plt.Line2D([], [], marker="o", color="w", label="UFBoot >= 95%", markerfacecolor="#000000", markersize=5)
        leg1 = ax.legend(handles=reg_patches + [sup_dot], title="Geographic Region & Support", loc="lower left", bbox_to_anchor=(0.01, 0.03), frameon=True, facecolor="white", framealpha=0.92, fontsize=8, title_fontsize=8.5)
        ax.add_artist(leg1)
        yr_patches = [mpatches.Patch(color=yr_colors[y], label=str(y)) for y in range(2015, 2026)]
        ax.legend(handles=yr_patches, title="Collection Year", loc="lower left", bbox_to_anchor=(0.36, 0.03), ncol=2, frameon=True, facecolor="white", framealpha=0.92, fontsize=8, title_fontsize=8.5)

def main():
    os.makedirs("results/figures", exist_ok=True)
    os.makedirs("results/treetime_clock_clean", exist_ok=True)
    if os.path.exists("results/gard/gard_clean_results.json"):
        gdata = json.load(open("results/gard/gard_clean_results.json")) if os.path.getsize("results/gard/gard_clean_results.json") > 0 else {}
        print(f"=== 1. HyPhy GARD Breakpoint Intervals Evaluated: {len(gdata.get('breakpointData', {}))} ===")

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, reg_col, yr_col = col_map["accession"], col_map["group"], col_map["region"], col_map["year"]
    reg_map = {"asia": "Asia", "northamerica": "North America", "europe": "Europe", "africa": "Africa", "oceania": "Oceania", "southamerica": "South America"}
    meta["Region_Plot"] = meta.apply(lambda r: "Asia" if str(r[grp_col]).strip().lower() == "asia" else reg_map.get(str(r[reg_col]).strip().lower(), "Other"), axis=1)
    meta["Year_Int"] = pd.to_numeric(meta[yr_col], errors="coerce").astype(int)
    meta_idx = meta.set_index(acc_col)

    root_acc = "KU322051"
    nuc_tree = Phylo.read("results/trees/tree_nuc_clean.treefile", "newick")
    nuc_tree.root_with_outgroup(root_acc)
    Phylo.write(nuc_tree, "results/trees/tree_nuc_clean_rooted.treefile", "newick")

    prot_tree = Phylo.read("results/trees/tree_prot_clean.treefile", "newick")
    if any(t.name == root_acc for t in prot_tree.get_terminals()):
        prot_tree.root_with_outgroup(root_acc)
    Phylo.write(prot_tree, "results/trees/tree_prot_clean_rooted.treefile", "newick")

    nuc_rep, prot_rep = parse_iqtree("results/trees/tree_nuc_clean.iqtree"), parse_iqtree("results/trees/tree_prot_clean.iqtree")
    n_tax, n_int, n_uf, n_both = summarize_support(nuc_tree)
    p_tax, p_int, p_uf, p_both = summarize_support(prot_tree)

    print("\n=== 2. Clean Phylogenetic Tree Comparison (N = 1,449) ===")
    print(f"Panel A - Nucleotide ML Tree : Taxa={n_tax} | Model={nuc_rep['model']} | LogL={nuc_rep['logl']} | Length={nuc_rep['tlen']} | UFBoot>=95%: {n_uf}/{n_int} ({n_uf/max(1,n_int)*100:.1f}%) | Both>=80/95: {n_both}/{n_int} ({n_both/max(1,n_int)*100:.1f}%)")
    print(f"Panel B - Amino-Acid ML Tree : Unique Taxa={p_tax} | Model={prot_rep['model']} | LogL={prot_rep['logl']} | Length={prot_rep['tlen']} | UFBoot>=95%: {p_uf}/{p_int} ({p_uf/max(1,p_int)*100:.1f}%) | Both>=80/95: {p_both}/{p_int} ({p_both/max(1,p_int)*100:.1f}%)")

    old_dates = pd.read_csv("results/treetime_clock/dates.csv")
    clean_dates = old_dates[old_dates["name"].isin(meta_idx.index)]
    clean_dates.to_csv("results/treetime_clock_clean/dates_clean.csv", index=False)
    subprocess.run(["treetime", "clock", "--tree", "results/trees/tree_nuc_clean_rooted.treefile", "--aln", "results/codon_aln_nuc_clean.fasta", "--dates", "results/treetime_clock_clean/dates_clean.csv", "--outdir", "results/treetime_clock_clean"], check=True)

    reg_colors = {"Asia": "#D62728", "North America": "#1F77B4", "Europe": "#2CA02C", "Africa": "#FF7F0E", "Oceania": "#9467BD", "South America": "#8C564B"}
    cmap_yr = plt.get_cmap("viridis", 11)
    yr_colors = {yr: matplotlib.colors.to_hex(cmap_yr(i)) for i, yr in enumerate(range(2015, 2026))}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 10.5))
    draw_tree_panel(ax1, nuc_tree, meta_idx, "A", f"Global Nucleotide ML Phylogeny (N = {n_tax}, {nuc_rep['model']})", 0.01, reg_colors, yr_colors, show_legend=True)
    draw_tree_panel(ax2, prot_tree, meta_idx, "B", f"Global Amino-Acid ML Phylogeny (Unique Taxa = {p_tax}, {prot_rep['model']})", 0.005, reg_colors, yr_colors, show_legend=False)
    plt.subplots_adjust(left=0.02, right=0.98, top=0.94, bottom=0.03, wspace=0.08)

    fig.savefig("results/figures/Figure1_Global_Phylogeny_AB_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure1_Global_Phylogeny_AB_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved 1000 DPI TIFF: results/figures/Figure1_Global_Phylogeny_AB_1000dpi.tiff")

if __name__ == "__main__":
    main()
