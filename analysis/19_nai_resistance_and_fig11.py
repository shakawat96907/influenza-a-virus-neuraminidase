#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO, Phylo

# Canonical WHO/CDC N1 Neuraminidase Inhibitor (NAI) Resistance Panel (N1 numbering)
NAI_PANEL = [
    ("E119V", 119, "E", "V", "Reduced inhibition (Oseltamivir / Peramivir)"),
    ("E119G", 119, "E", "G", "Reduced inhibition (Zanamivir / Peramivir)"),
    ("E119D", 119, "E", "D", "Reduced inhibition (Zanamivir / Laninamivir)"),
    ("Q136K", 136, "Q", "K", "Reduced inhibition (Zanamivir / Laninamivir)"),
    ("Y155H", 155, "Y", "H", "Reduced inhibition (Oseltamivir / Zanamivir)"),
    ("D199G", 199, "D", "G", "Reduced inhibition (Oseltamivir)"),
    ("D199E", 199, "D", "E", "Reduced inhibition (Oseltamivir)"),
    ("I223R", 223, "I", "R", "Reduced inhibition (Oseltamivir / Zanamivir; dual synergy with H275Y)"),
    ("I223V", 223, "I", "V", "Permissive / compensatory substitution (Oseltamivir)"),
    ("I223K", 223, "I", "K", "Reduced inhibition (Oseltamivir)"),
    ("I223T", 223, "I", "T", "Permissive / low-level reduced inhibition"),
    ("S247N", 247, "S", "N", "Reduced inhibition (Oseltamivir; dual synergy with H275Y)"),
    ("S247G", 247, "S", "G", "Reduced inhibition (Oseltamivir)"),
    ("H275Y", 275, "H", "Y", "Highly reduced inhibition (Oseltamivir / Peramivir)"),
    ("R293K", 293, "R", "K", "Highly reduced inhibition (Zanamivir / Oseltamivir)"),
    ("N295S", 295, "N", "S", "Reduced inhibition (Oseltamivir)")
]

def audit_phylogenetic_clustering(tree_path, target_accs):
    if len(target_accs) <= 1 or not os.path.exists(tree_path):
        return {acc: "Singleton" for acc in target_accs}, len(target_accs)
    tree = Phylo.read(tree_path, "newick")
    target_set = set(target_accs)
    monophyletic_groups = list()
    assigned = set()
    for clade in tree.find_clades(order="postorder"):
        tips = [t.name for t in clade.get_terminals()]
        if len(tips) >= 2 and set(tips).issubset(target_set):
            monophyletic_groups = [g for g in monophyletic_groups if not g.issubset(set(tips))]
            monophyletic_groups.append(set(tips))
    status_map = {}
    n_independent = 0
    for grp in monophyletic_groups:
        n_independent += 1
        for acc in grp:
            status_map[acc] = f"Monophyletic cluster (n={len(grp)})"
            assigned.add(acc)
    for acc in target_accs:
        if acc not in assigned:
            n_independent += 1
            status_map[acc] = "Independent singleton branch"
    return status_map, n_independent

def main():
    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)

    n1_to_aln = {}
    n1 = 0
    for idx, ch in enumerate(cons_gapped):
        if ch != "-":
            n1 += 1
            n1_to_aln[n1] = idx

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, yr_col = col_map["accession"], col_map["group"], col_map["year"]
    cnt_col = col_map.get("country", None)
    meta = meta[meta[acc_col].isin(aln_recs)].copy()
    meta["Is_Asia"] = meta[grp_col].str.strip().str.lower() == "asia"
    meta["Year_Int"] = pd.to_numeric(meta[yr_col], errors="coerce").astype(int)

    n_total = len(meta)
    n_asia = int(meta["Is_Asia"].sum())
    n_non = int((~meta["Is_Asia"]).sum())

    summary_rows = list()
    isolate_rows = list()

    for mut_label, pos, wt_aa, mut_aa, pheno in NAI_PANEL:
        c_idx = n1_to_aln[pos]
        hits_all, hits_asia, hits_non = list(), list(), list()

        for _, r in meta.iterrows():
            acc = r[acc_col]
            aa = aln_recs[acc][c_idx]
            if aa == mut_aa:
                hits_all.append(acc)
                if r["Is_Asia"]:
                    hits_asia.append(acc)
                else:
                    hits_non.append(acc)
                isolate_rows.append({
                    "Substitution": mut_label,
                    "N1_Position": pos,
                    "Accession": acc,
                    "Region": "Asia" if r["Is_Asia"] else "Non-Asia",
                    "Country": r[cnt_col] if cnt_col else "Unknown",
                    "Year": int(r["Year_Int"])
                })

        det_status = "Detected" if len(hits_all) > 0 else "Screened (0 detected)"
        summary_rows.append({
            "Substitution": mut_label,
            "N1_Position": pos,
            "Screening_Status": det_status,
            "Global_Count": len(hits_all),
            "Global_Pct": round(len(hits_all) / n_total * 100, 3),
            "Asian_Count": len(hits_asia),
            "Asian_Pct": round(len(hits_asia) / n_asia * 100, 3),
            "NonAsian_Count": len(hits_non),
            "NonAsian_Pct": round(len(hits_non) / n_non * 100, 3),
            "WHO_Phenotypic_Note": pheno
        })

    t14 = pd.DataFrame(summary_rows)
    t14.to_csv("results/tables/Main_Table_14_NAI_Resistance_Screening.csv", index=False)

    print("=================================================================")
    print(" MAIN TABLE 14: Screened vs. Detected NAI Resistance Substitutions")
    print("=================================================================")
    print(t14[["Substitution", "N1_Position", "Screening_Status", "Global_Count", "Global_Pct", "Asian_Count", "Asian_Pct", "NonAsian_Count", "NonAsian_Pct"]].to_string(index=False))

    iso_df = pd.DataFrame(isolate_rows)
    if len(iso_df) > 0:
        phylo_notes = list()
        for _, r in iso_df.iterrows():
            same_mut_accs = iso_df[iso_df["Substitution"] == r["Substitution"]]["Accession"].tolist()
            smap, _ = audit_phylogenetic_clustering("results/trees/tree_nuc_clean_rooted.treefile", same_mut_accs)
            phylo_notes.append(smap.get(r["Accession"], "Singleton"))
        iso_df["Phylogenetic_Clustering"] = phylo_notes
        iso_df = iso_df.sort_values(["N1_Position", "Substitution", "Region", "Year", "Country"])
        iso_df.to_csv("results/tables/Supplementary Table 16 - (detected_NAI_resistance_isolates_audit).csv", index=False)

        print("\n=================================================================")
        print(" SUPPLEMENTARY TABLE 16: Isolate-Level Audit of Detected Variants")
        print("=================================================================")
        print(iso_df.to_string(index=False))

        for sub_name in iso_df["Substitution"].unique():
            sub_df = iso_df[iso_df["Substitution"] == sub_name]
            _, n_ind = audit_phylogenetic_clustering("results/trees/tree_nuc_clean_rooted.treefile", sub_df["Accession"].tolist())
            print(f"\n[{sub_name}] Total={len(sub_df)} | Asia={(sub_df['Region']=='Asia').sum()} | Non-Asia={(sub_df['Region']=='Non-Asia').sum()} | Independent phylogenetic origins={n_ind}")

    # Render Clean, Overlap-Free Figure 11 (Legends Outside Axes; 300 DPI -> 1000 DPI metadata)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.2, 6.2), gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.25})

    x = np.arange(len(t14))
    w = 0.38
    b1 = ax1.bar(x - w/2, t14["Asian_Pct"], width=w, color="#C0392B", label="Asian Cohort (n = 750)")
    b2 = ax1.bar(x + w/2, t14["NonAsian_Pct"], width=w, color="#1F618D", label="Non-Asian Control (n = 699)")

    for idx_i, (_, r) in enumerate(t14.iterrows()):
        max_h = max(r["Asian_Pct"], r["NonAsian_Pct"])
        if r["Global_Count"] > 0:
            ax1.text(idx_i, max_h + 0.08, f"n={int(r['Global_Count'])}", ha="center", va="bottom", fontsize=7.5, fontweight="bold", color="#111111")
        else:
            ax1.text(idx_i, 0.05, "0", ha="center", va="bottom", fontsize=7.0, color="#7F8C8D")

    ax1.set_xticks(x)
    ax1.set_xticklabels(t14["Substitution"], rotation=45, ha="right", fontsize=8.8, fontweight="bold")
    ax1.set_ylabel("Prevalence (%)", fontsize=10, fontweight="bold")
    ax1.set_title("A   Screened vs. Detected Canonical N1 NAI Resistance Markers", loc="left", fontsize=10.5, fontweight="bold", pad=46)
    ax1.set_ylim(0, max(t14["Asian_Pct"].max(), t14["NonAsian_Pct"].max()) * 1.25 + 0.3)
    ax1.legend(handles=[b1, b2], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=8.2, frameon=False)
    ax1.grid(axis="y", linestyle=":", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    years = list(range(2015, 2026))
    det_subs = t14[t14["Global_Count"] > 0]["Substitution"].tolist()
    palette = ["#C0392B", "#2E86C1", "#27AE60", "#F39C12", "#8E44AD", "#16A085"]
    bottoms = np.zeros(len(years))
    handles_b = list()

    for s_idx, sub in enumerate(det_subs):
        counts_yr = [int(((iso_df["Substitution"] == sub) & (iso_df["Year"] == y)).sum()) for y in years]
        bars = ax2.bar(years, counts_yr, bottom=bottoms, width=0.62, color=palette[s_idx % len(palette)], label=f"{sub} (n={sum(counts_yr)})")
        handles_b.append(bars)
        bottoms += np.array(counts_yr)

    for idx_y, tot_y in enumerate(bottoms):
        if tot_y > 0:
            ax2.text(years[idx_y], tot_y + 0.15, f"{int(tot_y)}", ha="center", va="bottom", fontsize=8.2, fontweight="bold")

    ax2.set_xticks(years)
    ax2.set_xticklabels([str(y) for y in years], rotation=45, fontsize=8.8, fontweight="bold")
    ax2.set_xlabel("Collection Year (2015–2025)", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Number of Detected Isolates", fontsize=10, fontweight="bold")
    ax2.set_title("B   Temporal Distribution of Detected Resistance Substitutions", loc="left", fontsize=10.5, fontweight="bold", pad=46)

    ax2.set_ylim(0, max(bottoms.max() * 1.25, 4))
    ax2.legend(handles=handles_b, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=7.8, frameon=False)
    ax2.grid(axis="y", linestyle=":", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.055, right=0.98, top=0.77, bottom=0.16)
    png_path = "results/figures/Figure 11/Figure11_NAI_Resistance_Audit_1000dpi.png"
    tiff_path = "results/figures/Figure 11/Figure11_NAI_Resistance_Audit_1000dpi.tiff"
    fig.savefig(png_path, dpi=300, format="png")
    fig.savefig(tiff_path, dpi=300, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    for ext, pth in [("png", png_path), ("tiff", tiff_path)]:
        im = Image.open(pth)
        if ext == "tiff":
            im.save(pth, dpi=(1000, 1000), compression="tiff_lzw")
        else:
            im.save(pth, dpi=(1000, 1000))

    print("\nSaved New Figure 11 (1000 DPI TIFF & 1000 DPI PNG) in results/figures/Figure 11/ successfully!")

if __name__ == "__main__":
    main()
