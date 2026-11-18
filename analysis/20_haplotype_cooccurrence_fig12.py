#!/usr/bin/env python3
import itertools
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

VARIANT_PANEL = [
    ("M19T", 19, "T"), ("Q51K", 51, "K"), ("Y66F", 66, "F"), ("S74F", 74, "F"),
    ("R77G", 77, "G"), ("A81V", 81, "V"), ("S200N", 200, "N"), ("N222K", 222, "K"),
    ("I264T", 264, "T"), ("I389K", 389, "K"), ("N416D", 416, "D"), ("T452I", 452, "I"), ("V453M", 453, "M")
]
ASIA_ENRICHED = ["M19T", "Y66F", "R77G", "A81V", "S200N", "N222K", "V453M"]

def greedy_cover(series, target=0.90):
    tot = len(series)
    vc = series.value_counts()
    cum = 0
    rows = []
    for r, (hap, cnt) in enumerate(vc.items(), 1):
        cum += cnt
        rows.append({"Rank": r, "Haplotype": hap, "Asian_Count": int(cnt), "Asian_Pct": round(cnt/tot*100, 2), "Cum_Count": int(cum), "Cum_Pct": round(cum/tot*100, 2)})
        if cum / tot >= target and r >= 10:
            break
    n90 = next((x["Rank"] for x in rows if x["Cum_Pct"] >= target * 100), len(vc))
    return pd.DataFrame(rows), n90, len(vc)

def main():
    aln = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    n1_map = {}
    pos = 0
    for i, ch in enumerate(cons):
        if ch != "-":
            pos += 1
            n1_map[pos] = i

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    cm = {c.lower(): c for c in meta.columns}
    meta = meta[meta[cm["accession"]].isin(aln)].copy()
    meta["Is_Asia"] = meta[cm["group"]].str.strip().str.lower() == "asia"
    meta["Yr"] = pd.to_numeric(meta[cm["year"]], errors="coerce").astype(int)

    df = pd.DataFrame({"Acc": meta[cm["accession"]].values, "Is_Asia": meta["Is_Asia"].values, "Yr": meta["Yr"].values})
    for lbl, p, aa in VARIANT_PANEL:
        df[lbl] = [1 if aln[a][n1_map[p]] == aa else 0 for a in df["Acc"]]

    all_lbls = [v[0] for v in VARIANT_PANEL]
    df["Hap7"] = df.apply(lambda r: "+".join([c for c in ASIA_ENRICHED if r[c] == 1]) or "Consensus_Backbone", axis=1)
    df["Hap13"] = df.apply(lambda r: "+".join([c for c in all_lbls if r[c] == 1]) or "Consensus_Backbone", axis=1)

    asia = df[df["Is_Asia"]].copy()
    non = df[~df["Is_Asia"]].copy()
    print(f"All 7 Asia-enriched in 1 strain: {(asia[ASIA_ENRICHED].sum(axis=1)==7).sum()}/750 | Max in 1 strain: {asia[ASIA_ENRICHED].sum(axis=1).max()}/7")

    c7, n90_7, u7 = greedy_cover(asia["Hap7"], 0.90)
    c13, n90_13, u13 = greedy_cover(asia["Hap13"], 0.90)
    nc7 = non["Hap7"].value_counts()
    c7["NonAsian_Count"] = [int(nc7.get(h, 0)) for h in c7["Haplotype"]]
    c7["NonAsian_Pct"] = [round(nc7.get(h, 0)/len(non)*100, 2) for h in c7["Haplotype"]]
    c7["Years_Asia"] = [f"{asia[asia['Hap7']==h]['Yr'].min()}-{asia[asia['Hap7']==h]['Yr'].max()}" for h in c7["Haplotype"]]
    c7.to_csv("results/tables/Main_Table_15_Asian_Haplotype_SetCover.csv", index=False)

    print(f"7-Var Panel: Unique={u7}, Top1={c7['Asian_Pct'].iloc[0]}%, >=90% needs {n90_7} strains")
    print(f"13-Var Panel: Unique={u13}, Top1={c13['Asian_Pct'].iloc[0]}%, >=90% needs {n90_13} strains")
    print(c7.head(10).to_string(index=False))

    pairs, pvals = [], []
    for v1, v2 in itertools.combinations(all_lbls, 2):
        a = int(((asia[v1]==1) & (asia[v2]==1)).sum())
        b = int(((asia[v1]==1) & (asia[v2]==0)).sum())
        c = int(((asia[v1]==0) & (asia[v2]==1)).sum())
        d = int(((asia[v1]==0) & (asia[v2]==0)).sum())
        den = np.sqrt((a+b)*(c+d)*(a+c)*(b+d))
        phi = (a*d - b*c)/den if den > 0 else 0.0
        jacc = a/(a+b+c) if (a+b+c) > 0 else 0.0
        odds, pv = fisher_exact([[a, b], [c, d]])
        pvals.append(pv)
        pairs.append({"V1": v1, "V2": v2, "Both_Asia": a, "Phi": round(phi, 4), "Jaccard": round(jacc, 4), "OR": round(odds, 3) if np.isfinite(odds) else 999.0, "P": pv})

    pdf = pd.DataFrame(pairs)
    pdf["FDR"] = multipletests(pvals, method="fdr_bh")[1]
    pdf = pdf.sort_values("Phi", ascending=False)
    pdf.to_csv("results/tables/Supplementary Table 17 - (variant_pairwise_cooccurrence_network).csv", index=False)
    print("\nTop Positive Pairs:\n", pdf.head(6)[["V1", "V2", "Both_Asia", "Phi", "Jaccard", "FDR"]].to_string(index=False))
    print("\nTop Exclusive Pairs:\n", pdf.tail(5)[["V1", "V2", "Both_Asia", "Phi", "Jaccard", "FDR"]].to_string(index=False))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6.2), gridspec_kw={"width_ratios": [1.05, 1.15], "wspace": 0.30})
    nv = len(all_lbls)
    pm = np.eye(nv)
    im_idx = {l: i for i, l in enumerate(all_lbls)}
    for _, r in pdf.iterrows():
        i, j = im_idx[r["V1"]], im_idx[r["V2"]]
        pm[i, j] = pm[j, i] = r["Phi"]
    hm = ax1.imshow(pm, cmap="RdBu_r", vmin=-0.65, vmax=1.0, aspect="auto")
    ax1.set_xticks(range(nv)); ax1.set_yticks(range(nv))
    ax1.set_xticklabels(all_lbls, rotation=45, ha="right", fontsize=8.5, fontweight="bold")
    ax1.set_yticklabels(all_lbls, fontsize=8.5, fontweight="bold")
    ax1.set_title("A   Pairwise Variant Co-Occurrence Matrix (Phi φ, Asia n = 750)", loc="left", fontsize=10.2, fontweight="bold", pad=22)
    for i in range(nv):
        for j in range(nv):
            if i != j and abs(pm[i, j]) >= 0.25:
                ax1.text(j, i, f"{pm[i, j]:.2f}", ha="center", va="center", fontsize=6.5, fontweight="bold", color="white" if abs(pm[i, j]) > 0.55 else "black")
    fig.colorbar(hm, ax=ax1, fraction=0.046, pad=0.03).set_label("Phi Correlation (φ)", fontsize=8.8, fontweight="bold")

    kmax = 15
    v7 = [c7["Cum_Pct"].iloc[min(k, len(c7)-1)] for k in range(kmax)]
    v13 = [c13["Cum_Pct"].iloc[min(k, len(c13)-1)] for k in range(kmax)]
    rk = np.arange(1, kmax + 1)
    l1, = ax2.plot(rk, v7, marker="o", lw=2.2, color="#C0392B", label=f"Asia-Enriched 7-Variant Panel (≥90% needs {n90_7} strains)")
    l2, = ax2.plot(rk, v13, marker="s", lw=2.0, ls="--", color="#1F618D", label=f"Full 13-Variant Panel (≥90% needs {n90_13} strains)")
    ax2.axhline(90.0, color="#1E8449", ls=":", lw=1.2)
    ax2.set_xlabel("Distinct Haplotype Strains Included (Greedy Set Cover)", fontsize=9.8, fontweight="bold")
    ax2.set_ylabel("Cumulative Asian Cohort Coverage (%)", fontsize=9.8, fontweight="bold")
    ax2.set_title("B   Multi-Strain Haplotype Coverage Across Asian Isolates (n = 750)", loc="left", fontsize=10.2, fontweight="bold", pad=40)
    ax2.set_xticks(rk); ax2.set_ylim(0, 105)
    ax2.legend(handles=[l1, l2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=1, fontsize=8.0, frameon=False)
    ax2.grid(True, ls=":", alpha=0.35); ax2.spines["top"].set_visible(False); ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.98, top=0.79, bottom=0.14)
    p_png = "results/figures/Figure 12/Figure12_Haplotype_Cooccurrence_SetCover_1000dpi.png"
    p_tif = "results/figures/Figure 12/Figure12_Haplotype_Cooccurrence_SetCover_1000dpi.tiff"
    fig.savefig(p_png, dpi=300); fig.savefig(p_tif, dpi=300, pil_kwargs={"compression": "tiff_lzw"}); plt.close(fig)
    for ext, pth in [("png", p_png), ("tiff", p_tif)]:
        im = Image.open(pth)
        im.save(pth, dpi=(1000, 1000), compression="tiff_lzw") if ext == "tiff" else im.save(pth, dpi=(1000, 1000))
    print("SUCCESS: Saved Figure 12 (1000 DPI TIFF & PNG)!")

if __name__ == "__main__":
    main()
