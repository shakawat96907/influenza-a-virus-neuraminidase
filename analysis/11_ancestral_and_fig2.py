#!/usr/bin/env python3
import os, subprocess
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo, SeqIO
from Bio.Seq import Seq

def build_alignment_to_n1_map(cons_gapped):
    aln_to_n1 = {}
    n1 = 0
    for idx, ch in enumerate(cons_gapped):
        if ch != "-":
            n1 += 1
            aln_to_n1[idx] = (n1, ch)
    return aln_to_n1

def main():
    os.makedirs("results/ancestral", exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/figures", exist_ok=True)

    print("[1/3] Running TreeTime ancestral reconstruction (N = 1,449)...")
    subprocess.run([
        "treetime", "ancestral",
        "--tree", "results/trees/tree_nuc_clean_rooted.treefile",
        "--aln", "results/codon_aln_nuc_clean.fasta",
        "--outdir", "results/ancestral"
    ], check=True)

    print("[2/3] Counting parent -> child amino-acid transitions per site...")
    tree = Phylo.read("results/ancestral/annotated_tree.nexus", "nexus")
    anc_seqs = {r.id: str(r.seq) for r in SeqIO.parse("results/ancestral/ancestral_sequences.fasta", "fasta")}
    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    aln_to_n1 = build_alignment_to_n1_map(cons_gapped)

    prot_seqs = {}
    for node_id, nuc in anc_seqs.items():
        aa = []
        for i in range(0, len(nuc), 3):
            c = nuc[i:i+3]
            aa.append("-" if ("-" in c or len(c) < 3) else str(Seq(c).translate()))
        prot_seqs[node_id] = "".join(aa)

    site_events, sub_events = Counter(), Counter()
    sub_max_clade = defaultdict(int)
    tip_counts = {c: len(c.get_terminals()) for c in tree.find_clades()}

    for parent in tree.get_nonterminals():
        if parent.name not in prot_seqs: continue
        p_seq = prot_seqs[parent.name]
        for child in parent.clades:
            if child.name not in prot_seqs: continue
            c_seq = prot_seqs[child.name]
            n_desc = tip_counts[child]
            for aln_idx, (n1_pos, _) in aln_to_n1.items():
                pa, ca = p_seq[aln_idx], c_seq[aln_idx]
                if pa != "-" and ca != "-" and pa != ca and ca != "X":
                    site_events[n1_pos] += 1
                    sk = f"{pa}{n1_pos}{ca}"
                    sub_events[sk] += 1
                    if n_desc > sub_max_clade[sk]: sub_max_clade[sk] = n_desc

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col = col_map["accession"], col_map["group"]
    asia_ids = set(meta[meta[grp_col].str.lower() == "asia"][acc_col])
    nonasia_ids = set(meta[meta[grp_col].str.lower() != "asia"][acc_col])
    asia_ids, nonasia_ids = asia_ids, nonasia_ids
    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}

    rows = []
    for aln_idx, (n1_pos, cons_aa) in sorted(aln_to_n1.items(), key=lambda x: x):
        a_var = sum(1 for aid in (asia_ids & set(aln_recs)) if aln_recs[aid][aln_idx] not in ("-", cons_aa))
        na_var = sum(1 for nid in (nonasia_ids & set(aln_recs)) if aln_recs[nid][aln_idx] not in ("-", cons_aa))
        tot_var = a_var + na_var
        ind_orig = site_events[n1_pos]

        vc = Counter(aln_recs[sid][aln_idx] for sid in aln_recs if aln_recs[sid][aln_idx] not in ("-", cons_aa))
        dom_aa, _ = vc.most_common(1)[0] if vc else ("None", 0)
        dom_sub = f"{cons_aa}{n1_pos}{dom_aa}" if dom_aa != "None" else "None"
        dom_orig = sub_events.get(dom_sub, 0)
        max_clade = sub_max_clade.get(dom_sub, 0)

        if tot_var == 0: mech = "Invariant"
        elif dom_orig == 1 and max_clade >= 20: mech = "Single-origin founder"
        elif dom_orig >= 5 and max_clade >= 50: mech = "Trunk-sweeping & Recurrent"
        elif dom_orig >= 3: mech = "Recurrent / Homoplastic"
        else: mech = "Low-frequency / Sporadic"

        rows.append({
            "N1_Position": n1_pos, "Consensus_AA": cons_aa, "Dominant_Variant": dom_sub,
            "Global_Variant_Count": tot_var, "Global_Variant_Pct": round(tot_var / len(aln_recs) * 100, 2),
            "Asian_Variant_Pct": round(a_var / len(asia_ids) * 100, 2),
            "NonAsian_Variant_Pct": round(na_var / len(nonasia_ids) * 100, 2),
            "Site_Independent_Origins": ind_orig, "Dominant_Sub_Origins": dom_orig,
            "Max_Descendant_Clade_Size": max_clade, "Evolutionary_Mode": mech
        })

    df = pd.DataFrame(rows)
    df.to_csv("results/tables/Supplementary Table 8 - (ancestral_reconstruction_and_homoplasy).csv", index=False)

    print("\n=== Top 15 Most Frequent Variants & Independent Origins ===")
    print(df.sort_values("Global_Variant_Count", ascending=False).head(15)[["N1_Position", "Dominant_Variant", "Global_Variant_Pct", "Asian_Variant_Pct", "NonAsian_Variant_Pct", "Dominant_Sub_Origins", "Max_Descendant_Clade_Size", "Evolutionary_Mode"]].to_string(index=False))

    print("\n=== Top 15 Most Recurrent Sites (Independent Branch Events) ===")
    print(df.sort_values("Site_Independent_Origins", ascending=False).head(15)[["N1_Position", "Dominant_Variant", "Site_Independent_Origins", "Global_Variant_Pct", "Evolutionary_Mode"]].to_string(index=False))

    print("\n[3/3] Rendering New Figure 2 (1000 DPI TIFF & 300 DPI PNG)...")
    fig = plt.figure(figsize=(14, 10.5))
    axes = fig.subplots(3, 1, gridspec_kw={"height_ratios": [1.1, 1.1, 1.15], "hspace": 0.38}); ax_it = iter(axes)

    # Panel A: Independent Substitution Events per N1 site
    ax1 = next(ax_it)
    ax1.bar(df["N1_Position"], df["Site_Independent_Origins"], width=1.1, color="#2B5C8F", alpha=0.9)
    ax1.axvspan(1, 35, color="#E0E0E0", alpha=0.4, label="TM / Cytoplasmic (1–35)")
    ax1.axvspan(36, 82, color="#F5DEB3", alpha=0.35, label="Stalk (36–82)")
    for _, r in df.nlargest(8, "Site_Independent_Origins").iterrows():
        ax1.text(r["N1_Position"], r["Site_Independent_Origins"] + 0.8, f"{r['Dominant_Variant']} ({int(r['Site_Independent_Origins'])})", fontsize=7.5, ha="center", fontweight="bold")
    ax1.set_ylabel("Independent Events", fontsize=10, fontweight="bold")
    ax1.set_title("A   Phylogenetically Independent Amino-Acid Substitution Events per Site (TreeTime Ancestral)", loc="left", fontsize=11, fontweight="bold")
    ax1.set_xlim(0, 472)
    ax1.set_ylim(0, df["Site_Independent_Origins"].max() * 1.18)
    ax1.legend(loc="upper right", fontsize=8, frameon=True)
    ax1.grid(axis="y", linestyle="--", alpha=0.3)

    # Panel B: Non-Consensus Variant Frequency (%) Asia vs Non-Asia
    ax2 = next(ax_it)
    ax2.plot(df["N1_Position"], df["Asian_Variant_Pct"], color="#D62728", lw=1.0, alpha=0.85, label="Asian Cohort (n = 750)")
    ax2.plot(df["N1_Position"], -df["NonAsian_Variant_Pct"], color="#1F77B4", lw=1.0, alpha=0.85, label="Non-Asian Control (n = 699)")
    ax2.axhline(0, color="black", lw=0.8)
    for _, r in df.nlargest(8, "Global_Variant_Pct").iterrows():
        ax2.text(r["N1_Position"], r["Asian_Variant_Pct"] + 3, r["Dominant_Variant"], fontsize=7.5, ha="center", color="#8B0000", fontweight="bold")
    ax2.set_ylabel("Variant Freq. (%)\n(Non-Asia | Asia)", fontsize=9.5, fontweight="bold")
    ax2.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Non-Consensus Amino-Acid Variant Prevalence Across Asian and Non-Asian Cohorts", loc="left", fontsize=11, fontweight="bold")
    ax2.set_xlim(0, 472)
    ax2.set_ylim(-115, 118)
    ax2.set_yticks([-100, -50, 0, 50, 100])
    ax2.set_yticklabels(["100%", "50%", "0%", "50%", "100%"])
    ax2.legend(loc="upper right", fontsize=8.5, frameon=True)
    ax2.grid(axis="y", linestyle="--", alpha=0.3)

    # Panel C: Top 18 Polymorphic Sites - Raw Frequency vs Independent Origins
    ax3 = next(ax_it)
    top18 = df.nlargest(18, "Global_Variant_Count").copy()
    x_pos = np.arange(len(top18))
    ax3_twin = ax3.twinx()
    b1 = ax3.bar(x_pos - 0.18, top18["Global_Variant_Pct"], width=0.35, color="#E67E22", label="Global Variant Frequency (%)")
    b2 = ax3_twin.bar(x_pos + 0.18, top18["Site_Independent_Origins"], width=0.35, color="#2E86C1", label="Independent Origins (Events)")
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(top18["Dominant_Variant"], rotation=35, ha="right", fontsize=8.5, fontweight="bold")
    ax3.set_ylabel("Global Frequency (%)", fontsize=9.5, fontweight="bold", color="#D35400")
    ax3_twin.set_ylabel("Independent Events", fontsize=9.5, fontweight="bold", color="#1B4F72")
    ax3.set_title("C   Distinguishing Clade-Sweeping Variants from Recurrent Homoplastic Substitutions (Top 18 Sites)", loc="left", fontsize=11, fontweight="bold")
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc="upper right", fontsize=8.5, frameon=True)
    ax3.grid(axis="y", linestyle="--", alpha=0.3)

    fig.savefig("results/figures/Figure2_Ancestral_And_Variant_Landscape_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure2_Ancestral_And_Variant_Landscape_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("Saved New Figure 2 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
