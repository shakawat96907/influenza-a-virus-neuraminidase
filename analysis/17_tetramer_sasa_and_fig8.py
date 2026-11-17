#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, Structure, Model
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley

TIEN_MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0
}

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
}

def classify_rsa(rsa):
    if rsa >= 0.25:
        return "Highly surface-exposed (RSA ≥ 0.25)"
    elif rsa >= 0.05:
        return "Moderately exposed (0.05 ≤ RSA < 0.25)"
    return "Buried / Interfacial (RSA < 0.05)"

def main():
    pdb_path = "results/structure/4B7R_tetramer.pdb"
    parser = PDBParser(QUIET=True)
    raw_struct = parser.get_structure("4B7R_bio", pdb_path)

    tet_struct = Structure.Structure("tetramer")
    tet_model = Model.Model(0)
    tet_struct.add(tet_model)

    chain_letters = ["A", "B", "C", "D", "E", "F", "G", "H"]
    c_idx = 0
    for mdl in raw_struct:
        for ch in mdl:
            aa_res = [r for r in ch if is_aa(r, standard=True) and r.resname in TIEN_MAX_ASA]
            if len(aa_res) > 100:
                new_ch = ch.copy()
                new_ch.id = chain_letters[c_idx]
                for r in list(new_ch):
                    if not (is_aa(r, standard=True) and r.resname in TIEN_MAX_ASA):
                        new_ch.detach_child(r.get_id())
                tet_model.add(new_ch)
                c_idx += 1

    print(f"Loaded biological assembly with {c_idx} protein subunits. Computing SASA & RSA...")
    sr = ShrakeRupley()
    sr.compute(tet_struct, level="R")
    tet_sasa = {int(r.get_id()[1]): float(r.sasa) for r in tet_model["A"]}

    mono_struct = Structure.Structure("monomer")
    mono_model = Model.Model(0)
    mono_struct.add(mono_model)
    mono_model.add(tet_model["A"].copy())
    sr.compute(mono_struct, level="R")
    mono_sasa = {int(r.get_id()[1]): float(r.sasa) for r in mono_model["A"]}

    catalytic_set = {118, 119, 151, 152, 179, 180, 223, 225, 247, 275, 277, 293, 295, 368, 402}
    rows = []
    for r in mono_model["A"]:
        pos = int(r.get_id()[1])
        aa3 = r.resname
        aa1 = AA3_TO_1[aa3]
        max_asa = TIEN_MAX_ASA[aa3]
        m_sasa = mono_sasa[pos]
        t_sasa = tet_sasa[pos]
        d_sasa = m_sasa - t_sasa
        t_rsa = t_sasa / max_asa
        buried_int = "Yes (Interface)" if d_sasa >= 10.0 else "No"

        rows.append({
            "N1_Position": pos, "Residue_4B7R": f"{aa1}{pos}", "AA3": aa3,
            "MaxASA_Tien2013": max_asa,
            "Monomer_SASA_A2": round(m_sasa, 2),
            "Tetramer_SASA_A2": round(t_sasa, 2),
            "Delta_SASA_Buried": round(d_sasa, 2),
            "Monomer_RSA": round(m_sasa / max_asa, 4),
            "Tetramer_RSA": round(t_rsa, 4),
            "Interface_Buried": buried_int,
            "Tetramer_Exposure_Class": classify_rsa(t_rsa),
            "Role": "Catalytic / Active-Site" if pos in catalytic_set else "Head Domain"
        })

    df = pd.DataFrame(rows)
    df.to_csv("results/tables/Supplementary Table 14 - (tetramer_vs_monomer_SASA_RSA_4B7R).csv", index=False)

    n_int_buried = (df["Delta_SASA_Buried"] >= 10.0).sum()
    print(f"Total head residues (83–469): {len(df)} | Residues buried at tetramer interface (ΔSASA ≥ 10 Å²): {n_int_buried}")

    key_head = [118, 151, 152, 200, 222, 225, 264, 270, 275, 277, 293, 339, 368, 382, 389, 402, 416, 449, 452, 453, 469]
    t12 = df[df["N1_Position"].isin(key_head)].copy()
    t12.to_csv("results/tables/Main_Table_12_Tetramer_SASA_RSA_Candidates.csv", index=False)

    print("\n=================================================================")
    print(" MAIN TABLE 12: Monomer vs. Biological Tetramer SASA & RSA (4B7R)")
    print("=================================================================")
    print(t12[["N1_Position", "Residue_4B7R", "Role", "Monomer_SASA_A2", "Tetramer_SASA_A2", "Delta_SASA_Buried", "Tetramer_RSA", "Interface_Buried", "Tetramer_Exposure_Class"]].to_string(index=False))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6.8), gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.25})
    non_int = df[df["Delta_SASA_Buried"] < 10.0]
    yes_int = df[df["Delta_SASA_Buried"] >= 10.0]
    s1 = ax1.scatter(non_int["Monomer_SASA_A2"], non_int["Tetramer_SASA_A2"], s=18, color="#2E86C1", alpha=0.72, label=f"Non-interfacial (n = {len(non_int)})")
    s2 = ax1.scatter(yes_int["Monomer_SASA_A2"], yes_int["Tetramer_SASA_A2"], s=28, color="#C0392B", alpha=0.85, edgecolors="#111111", lw=0.4, label=f"Buried at Tetramer Interface (ΔSASA ≥ 10 Å², n = {len(yes_int)})")
    ax1.plot([0, 220], [0, 220], color="#333333", linestyle="--", lw=1.0)

    for _, r in yes_int.nlargest(5, "Delta_SASA_Buried").iterrows():
        ax1.annotate(r["Residue_4B7R"], xy=(r["Monomer_SASA_A2"], r["Tetramer_SASA_A2"]), xytext=(r["Monomer_SASA_A2"] + 8, max(2, r["Tetramer_SASA_A2"] - 4)), fontsize=8, fontweight="bold", color="#78281F", arrowprops=dict(arrowstyle="-", lw=0.6, color="#555555"))

    ax1.set_xlabel("Isolated Monomer SASA (Å², Chain A)", fontsize=10, fontweight="bold")
    ax1.set_ylabel("Biological Homotetramer SASA (Å², 4B7R)", fontsize=10, fontweight="bold")
    ax1.set_title("A   Monomer vs. Homotetramer Solvent Accessibility (PDB 4B7R)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax1.set_xlim(-5, 225)
    ax1.set_ylim(-5, 225)
    ax1.legend(handles=[s1, s2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8.0, frameon=False)
    ax1.grid(True, linestyle=":", alpha=0.4)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    exp_hi = df[df["Tetramer_RSA"] >= 0.25]
    exp_md = df[(df["Tetramer_RSA"] >= 0.05) & (df["Tetramer_RSA"] < 0.25)]
    exp_bu = df[df["Tetramer_RSA"] < 0.05]

    b1 = ax2.scatter(exp_bu["N1_Position"], exp_bu["Tetramer_RSA"], s=14, color="#95A5A6", alpha=0.65, label=f"Buried (RSA < 0.05, n = {len(exp_bu)})")
    b2 = ax2.scatter(exp_md["N1_Position"], exp_md["Tetramer_RSA"], s=20, color="#F39C12", alpha=0.80, label=f"Moderate (0.05–0.25, n = {len(exp_md)})")
    b3 = ax2.scatter(exp_hi["N1_Position"], exp_hi["Tetramer_RSA"], s=28, color="#27AE60", alpha=0.88, label=f"Exposed (≥ 0.25, n = {len(exp_hi)})")
    ax2.axhline(0.25, color="#1E8449", linestyle="--", lw=0.9)
    ax2.axhline(0.05, color="#7F8C8D", linestyle=":", lw=0.9)

    label_pos = [200, 222, 264, 339, 382, 389, 416, 452, 453, 469]
    for idx_k, (_, r) in enumerate(df[df["N1_Position"].isin(label_pos)].iterrows()):
        dy = 0.07 if idx_k % 2 == 0 else -0.09
        ax2.annotate(r["Residue_4B7R"], xy=(r["N1_Position"], r["Tetramer_RSA"]), xytext=(r["N1_Position"], min(0.95, max(0.03, r["Tetramer_RSA"] + dy))), fontsize=8, fontweight="bold", ha="center", color="#1A2530", arrowprops=dict(arrowstyle="-", lw=0.6, color="#555555"))

    ax2.set_xlabel("Canonical N1 Head-Domain Position (83–469)", fontsize=10, fontweight="bold")
    ax2.set_ylabel("Tetramer Relative Solvent Accessibility (RSA)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Biological Tetramer RSA Profile & Variant Surface Exposure", loc="left", fontsize=10.5, fontweight="bold", pad=26)
    ax2.set_xlim(78, 475)
    ax2.set_ylim(-0.04, 1.05)
    ax2.legend(handles=[b3, b2, b1], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=7.8, frameon=False)
    ax2.grid(True, linestyle=":", alpha=0.4)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    plt.subplots_adjust(left=0.06, right=0.97, top=0.86, bottom=0.12)
    fig.savefig("results/figures/Figure8_Tetramer_SASA_And_RSA_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure8_Tetramer_SASA_And_RSA_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("Saved New Figure 8 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
