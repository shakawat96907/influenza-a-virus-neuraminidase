#!/usr/bin/env python3
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO, Entrez

Entrez.email = "farjana-bmb@sust.edu"

def ensure_clean_who_vaccines():
    recs = list(SeqIO.parse("results/vaccine/who_vaccine_na_refs.fasta", "fasta"))
    audit = pd.read_csv("results/vaccine/vaccine_reference_audit.csv")

    # Replace 2020 if it points to non-pdm09 QFQ15359.1
    need_2020 = any(r.id.startswith("2020|QFQ15359") for r in recs)
    if need_2020:
        sh = Entrez.esearch(db="protein", term='txid11320[Organism] AND "Brisbane/02/2018"[Title] AND 469[SLEN]', retmax=5)
        ids = Entrez.read(sh).get("IdList", list())
        sh.close()
        time.sleep(0.4)
        if not ids:
            sh = Entrez.esearch(db="protein", term='txid11320[Organism] AND Brisbane[All Fields] AND 2018[All Fields] AND "MNPNQKIITIGSICMTIG"[All Fields] AND 469[SLEN]', retmax=5)
            ids = Entrez.read(sh).get("IdList", list())
            sh.close()
            time.sleep(0.4)
        fh = Entrez.efetch(db="protein", id=ids.pop(0), rettype="gb", retmode="text")
        gb20 = SeqIO.read(fh, "genbank")
        fh.close()
        gb20.id = f"2020|{gb20.id}"
        gb20.description = "A/Brisbane/02/2018 (2020 WHO Vaccine)"
        recs = [gb20 if r.id.startswith("2020|") else r for r in recs]
        SeqIO.write(recs, "results/vaccine/who_vaccine_na_refs.fasta", "fasta")

        audit.loc[audit["Season_Tag"] == "2020", "Verified_Accession"] = gb20.id.split("|").pop()
        audit.loc[audit["Season_Tag"] == "2020", "WHO_Label"] = "A/Brisbane/02/2018 (2020 WHO Vaccine)"
        audit.loc[audit["Season_Tag"] == "2020", "NCBI_Organism"] = gb20.annotations.get("organism", "Influenza A virus")
        audit.loc[audit["Season_Tag"] == "2020", "NCBI_Description"] = gb20.description
        audit.to_csv("results/vaccine/vaccine_reference_audit.csv", index=False)

    print("=== Final Verified WHO Vaccine & Legacy Reference Audit ===")
    print(audit[["Season_Tag", "Verified_Accession", "WHO_Label", "Length_AA"]].to_string(index=False))
    return {r.id.split("|").pop(0): str(r.seq) for r in recs}

def year_to_season_tag(y):
    if y <= 2016: return "2015_2016"
    elif y <= 2019: return "2017_2019"
    elif y == 2020: return "2020"
    elif y <= 2022: return "2021_2022"
    return "2023_2025"

def main():
    ref_seqs = ensure_clean_who_vaccines()
    # All 6 sequences are exact 469 aa canonical N1 length (1-to-1 positional index)
    for k, s in ref_seqs.items():
        assert len(s) == 469, f"Expected 469 aa for {k}, got {len(s)}"

    aln_recs = {r.id: str(r.seq) for r in SeqIO.parse("results/aln_prot_clean.fasta", "fasta")}
    cons_gapped = str(next(SeqIO.parse("results/consensus_prot_clean.fasta", "fasta")).seq)
    cons_nogap = cons_gapped.replace("-", "")

    n1_to_aln = {}
    n1 = 0
    for idx, ch in enumerate(cons_gapped):
        if ch != "-":
            n1 += 1
            n1_to_aln[n1] = idx

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, yr_col = col_map["accession"], col_map["group"], col_map["year"]
    meta = meta[meta[acc_col].isin(aln_recs)].copy()
    meta["Is_Asia"] = meta[grp_col].str.strip().str.lower() == "asia"
    meta["Year_Int"] = pd.to_numeric(meta[yr_col], errors="coerce").astype(int)
    meta["Season_Tag"] = meta["Year_Int"].apply(year_to_season_tag)

    asia_meta = meta[meta["Is_Asia"]]
    non_meta = meta[~meta["Is_Asia"]]

    rows = []
    for pos in range(1, 470):
        c_idx = n1_to_aln[pos]
        cons_aa = cons_nogap[pos - 1]
        legacy_aa = ref_seqs["2009_Legacy"][pos - 1]

        leg_div_all = sum(1 for a in meta[acc_col] if aln_recs[a][c_idx] not in ("-", "X", legacy_aa)) / len(meta) * 100
        leg_div_asia = sum(1 for a in asia_meta[acc_col] if aln_recs[a][c_idx] not in ("-", "X", legacy_aa)) / len(asia_meta) * 100

        sm_div_all = sum(
            1 for _, r in meta.iterrows()
            if aln_recs[r[acc_col]][c_idx] not in ("-", "X", ref_seqs[r["Season_Tag"]][pos - 1])
        ) / len(meta) * 100

        sm_div_asia = sum(
            1 for _, r in asia_meta.iterrows()
            if aln_recs[r[acc_col]][c_idx] not in ("-", "X", ref_seqs[r["Season_Tag"]][pos - 1])
        ) / len(asia_meta) * 100

        sm_div_non = sum(
            1 for _, r in non_meta.iterrows()
            if aln_recs[r[acc_col]][c_idx] not in ("-", "X", ref_seqs[r["Season_Tag"]][pos - 1])
        ) / len(non_meta) * 100

        rows.append({
            "N1_Position": pos, "Consensus_AA": cons_aa, "Legacy_QFQ15359_AA": legacy_aa,
            "Legacy_Ref_Div_Global_Pct": round(leg_div_all, 2),
            "Legacy_Ref_Div_Asia_Pct": round(leg_div_asia, 2),
            "SeasonMatched_Div_Global_Pct": round(sm_div_all, 2),
            "SeasonMatched_Div_Asia_Pct": round(sm_div_asia, 2),
            "SeasonMatched_Div_NonAsia_Pct": round(sm_div_non, 2),
            "Artifact_Reduction_Pct": round(leg_div_all - sm_div_all, 2)
        })

    df = pd.DataFrame(rows)
    df.to_csv("results/tables/Supplementary Table 15 - (season_matched_vaccine_divergence_469).csv", index=False)

    n_leg_gt50 = (df["Legacy_Ref_Div_Global_Pct"] >= 50.0).sum()
    n_sm_gt50 = (df["SeasonMatched_Div_Global_Pct"] >= 50.0).sum()
    print(f"\nSites with ≥50% divergence -> Legacy Ref (QFQ15359.1): {n_leg_gt50} sites | Season-Matched WHO Vaccines: {n_sm_gt50} sites")

    t13 = df[df["SeasonMatched_Div_Global_Pct"] >= 10.0].sort_values("SeasonMatched_Div_Global_Pct", ascending=False).copy()
    t13.to_csv("results/tables/Main_Table_13_SeasonMatched_Vaccine_Divergence.csv", index=False)

    print("\n=================================================================")
    print(" MAIN TABLE 13: Top Season-Matched WHO Vaccine Divergent Sites (≥10%)")
    print("=================================================================")
    print(t13.to_string(index=False))

    # Render Clean, Overlap-Free Figure 9 at 1000 DPI (both TIFF and PNG)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14.5, 9.0), gridspec_kw={"height_ratios": [1.05, 1.0], "hspace": 0.48})

    l1, = ax1.plot(df["N1_Position"], df["Legacy_Ref_Div_Global_Pct"], color="#95A5A6", lw=1.0, alpha=0.80, label=f"Misassigned Legacy Reference QFQ15359.1 ({n_leg_gt50} sites ≥50%)")
    l2, = ax1.plot(df["N1_Position"], df["SeasonMatched_Div_Global_Pct"], color="#C0392B", lw=1.45, alpha=0.95, label=f"Season-Matched WHO H1N1pdm09 Vaccines ({n_sm_gt50} sites ≥50%)")

    ax1.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Divergence from Vaccine Strain (%)", fontsize=10, fontweight="bold")
    ax1.set_title("A   Elimination of Non-pdm09 Reference Artifacts via Season-Matched WHO Vaccine Comparison (N = 1,449)", loc="left", fontsize=10.6, fontweight="bold", pad=30)
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
    ax2.set_title("B   Regional Comparison of Top Season-Matched WHO Vaccine-Divergent Sites (2015–2025)", loc="left", fontsize=10.6, fontweight="bold", pad=30)
    ax2.set_ylim(0, max(top_b["SeasonMatched_Div_Asia_Pct"].max(), top_b["SeasonMatched_Div_NonAsia_Pct"].max()) * 1.18)
    ax2.legend(handles=[b_a, b_na], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=8.4, frameon=False)
    ax2.grid(axis="y", linestyle=":", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.97, top=0.90, bottom=0.08)
    fig.savefig("results/figures/Figure9_SeasonMatched_Vaccine_Divergence_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure9_SeasonMatched_Vaccine_Divergence_1000dpi.png", dpi=1000, format="png")
    plt.close(fig)
    print("\nSaved New Figure 9 (1000 DPI TIFF & 1000 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
