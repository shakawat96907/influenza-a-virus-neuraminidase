#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    # 1. Generate Main Table 5 from T5 outputs
    t6 = pd.read_csv("results/tables/Supplementary Table 6 - (intercontinental_migration_transitions).csv")
    t7 = pd.read_csv("results/tables/Supplementary Table 7 - (asia_enriched_clades).csv")

    reg_names = {
        "North_America": "North America",
        "Europe": "Europe",
        "Oceania": "Oceania",
        "Africa": "Africa",
        "South_America": "South America"
    }
    table5_rows = []
    for reg_key, reg_label in reg_names.items():
        exp_cnt = int(t6[(t6["From_Region"] == "Asia") & (t6["To_Region"] == reg_key)]["Transition_Count"].sum())
        imp_cnt = int(t6[(t6["From_Region"] == reg_key) & (t6["To_Region"] == "Asia")]["Transition_Count"].sum())
        table5_rows.append({
            "Partner_Region": reg_label,
            "Exports_Asia_to_Region": exp_cnt,
            "Introductions_Region_to_Asia": imp_cnt,
            "Net_Flow_From_Asia": exp_cnt - imp_cnt,
            "Total_Bilateral_Transitions": exp_cnt + imp_cnt
        })
    t5_df = pd.DataFrame(table5_rows)
    t5_df.loc[len(t5_df)] = {
        "Partner_Region": "Total (Asia <-> Non-Asia)",
        "Exports_Asia_to_Region": t5_df["Exports_Asia_to_Region"].sum(),
        "Introductions_Region_to_Asia": t5_df["Introductions_Region_to_Asia"].sum(),
        "Net_Flow_From_Asia": t5_df["Net_Flow_From_Asia"].sum(),
        "Total_Bilateral_Transitions": t5_df["Total_Bilateral_Transitions"].sum()
    }
    t5_df.to_csv("results/tables/Main_Table_5_Migration_Summary.csv", index=False)

    print("=================================================================")
    print(" MAIN TABLE 5: Bidirectional Intercontinental Migration (N=1,449)")
    print("=================================================================")
    print(t5_df.to_string(index=False))
    print(f"\nNote: Plus 55 Non-Asia <-> Non-Asia transitions (Grand Total = 338) and {len(t7)} Asia-enriched subclusters (>=80% Asian; 32 with 100% Asian isolates).\n")

    # 2. Generate Main Table 6 from T6 output
    df = pd.read_csv("results/tables/Supplementary Table 8 - (ancestral_reconstruction_and_homoplasy).csv")
    top15 = df.sort_values("Global_Variant_Count", ascending=False).head(15).copy()
    t6_main = top15[[
        "N1_Position", "Dominant_Variant", "Global_Variant_Pct",
        "Asian_Variant_Pct", "NonAsian_Variant_Pct",
        "Site_Independent_Origins", "Max_Descendant_Clade_Size", "Evolutionary_Mode"
    ]].rename(columns={
        "N1_Position": "N1_Pos",
        "Dominant_Variant": "Variant",
        "Global_Variant_Pct": "Global_%",
        "Asian_Variant_Pct": "Asian_%",
        "NonAsian_Variant_Pct": "NonAsian_%",
        "Site_Independent_Origins": "Indep_Events",
        "Max_Descendant_Clade_Size": "Max_Clade",
        "Evolutionary_Mode": "Evolutionary_Mode"
    })
    t6_main.to_csv("results/tables/Main_Table_6_Top15_Variants_Ancestral.csv", index=False)

    print("=================================================================")
    print(" MAIN TABLE 6: Top 15 Variants & Independent Origins (TreeTime)")
    print("=================================================================")
    print(t6_main.to_string(index=False))

    # 3. Render Clean, Overlap-Free Figure 2 (Legends Outside Axes)
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(13.5, 11.5),
        gridspec_kw={"height_ratios": [1.05, 1.05, 1.15], "hspace": 0.48}
    )

    # Panel A: Independent Substitution Events per N1 site
    ax1.bar(df["N1_Position"], df["Site_Independent_Origins"], width=1.15, color="#2B5C8F", alpha=0.92, zorder=2)
    p_tm = ax1.axvspan(1, 35, color="#D9D9D9", alpha=0.45, label="TM / Cytoplasmic (1–35)", zorder=1)
    p_st = ax1.axvspan(36, 82, color="#F5DEB3", alpha=0.45, label="Stalk Region (36–82)", zorder=1)
    p_hd = ax1.axvspan(83, 469, color="#EBF3FA", alpha=0.35, label="Catalytic Head Domain (83–469)", zorder=1)

    # Staggered non-overlapping callouts for top recurrent sites
    top_a = df.nlargest(8, "Site_Independent_Origins").sort_values("N1_Position")
    offsets_a = [(14, 38), (10, 72), (18, 42), (-22, 45), (-28, 35), (26, 68), (-34, 42), (-18, 72)]
    for (_, r), (dx, dy) in zip(top_a.iterrows(), offsets_a):
        x, y = r["N1_Position"], r["Site_Independent_Origins"]
        ax1.annotate(
            f"{r['Dominant_Variant']} ({int(y)})",
            xy=(x, y), xytext=(x + dx, y + dy),
            fontsize=8, fontweight="bold", ha="center", color="#1A2530",
            arrowprops=dict(arrowstyle="-", color="#555555", lw=0.7, shrinkA=0, shrinkB=2)
        )

    ax1.set_ylabel("Independent Events", fontsize=10, fontweight="bold")
    ax1.set_title("A   Phylogenetically Independent Amino-Acid Substitution Events per Site (TreeTime Ancestral)", loc="left", fontsize=11, fontweight="bold", pad=24)
    ax1.set_xlim(-8, 485)
    ax1.set_ylim(0, df["Site_Independent_Origins"].max() * 1.28)
    ax1.legend(handles=[p_tm, p_st, p_hd], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=3, fontsize=8.5, frameon=False)
    ax1.grid(axis="y", linestyle="--", alpha=0.3, zorder=0)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    # Panel B: Non-Consensus Variant Frequency (%) Asia vs Non-Asia
    l_asia, = ax2.plot(df["N1_Position"], df["Asian_Variant_Pct"], color="#C0392B", lw=1.1, alpha=0.9, label="Asian Cohort (n = 750)")
    l_non, = ax2.plot(df["N1_Position"], -df["NonAsian_Variant_Pct"], color="#1F618D", lw=1.1, alpha=0.9, label="Non-Asian Control (n = 699)")
    ax2.axhline(0, color="black", lw=0.8)

    top_b = df.nlargest(8, "Global_Variant_Pct").sort_values("N1_Position")
    offsets_b = [(-18, 16), (0, 22), (16, 16), (-20, 18), (18, 22), (-16, 18), (18, 22), (14, 16)]
    for (_, r), (dx, dy) in zip(top_b.iterrows(), offsets_b):
        x, y = r["N1_Position"], r["Asian_Variant_Pct"]
        ax2.annotate(
            r["Dominant_Variant"],
            xy=(x, y), xytext=(x + dx, y + dy),
            fontsize=8, fontweight="bold", ha="center", color="#78281F",
            arrowprops=dict(arrowstyle="-", color="#922B21", lw=0.7, shrinkA=0, shrinkB=2)
        )

    ax2.set_ylabel("Variant Freq. (%)\n(Non-Asia | Asia)", fontsize=9.5, fontweight="bold")
    ax2.set_xlabel("Canonical N1 Amino-Acid Position (1–469)", fontsize=10, fontweight="bold")
    ax2.set_title("B   Non-Consensus Amino-Acid Variant Prevalence Across Asian and Non-Asian Cohorts", loc="left", fontsize=11, fontweight="bold", pad=24)
    ax2.set_xlim(0, 475)
    ax2.set_ylim(-110, 118)
    ax2.set_yticks([-100, -50, 0, 50, 100])
    ax2.set_yticklabels(["100%", "50%", "0%", "50%", "100%"])
    ax2.legend(handles=[l_asia, l_non], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=8.5, frameon=False)
    ax2.grid(axis="y", linestyle="--", alpha=0.3)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    # Panel C: Top 15 Polymorphic Sites - Raw Frequency vs Independent Origins
    x_pos = np.arange(len(top15))
    ax3_twin = ax3.twinx()
    b1 = ax3.bar(x_pos - 0.18, top15["Global_Variant_Pct"], width=0.35, color="#E67E22", label="Global Variant Frequency (%)")
    b2 = ax3_twin.bar(x_pos + 0.18, top15["Site_Independent_Origins"], width=0.35, color="#2E86C1", label="Independent Origins (Events)")
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(top15["Dominant_Variant"], rotation=32, ha="right", fontsize=8.5, fontweight="bold")
    ax3.set_ylabel("Global Frequency (%)", fontsize=9.5, fontweight="bold", color="#BA4A00")
    ax3_twin.set_ylabel("Independent Events", fontsize=9.5, fontweight="bold", color="#1B4F72")
    ax3.set_title("C   Clade-Sweeping Variants vs. Recurrent Homoplastic Substitutions (Top 15 Sites)", loc="left", fontsize=11, fontweight="bold", pad=24)
    ax3.set_ylim(0, top15["Global_Variant_Pct"].max() * 1.18)
    ax3_twin.set_ylim(0, top15["Site_Independent_Origins"].max() * 1.18)
    ax3.legend(handles=[b1, b2], loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=2, fontsize=8.5, frameon=False)
    ax3.grid(axis="y", linestyle="--", alpha=0.3)
    ax3.spines["top"].set_visible(False)
    ax3_twin.spines["top"].set_visible(False)

    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07)
    fig.savefig("results/figures/Figure2_Ancestral_And_Variant_Landscape_v2_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig("results/figures/Figure2_Ancestral_And_Variant_Landscape_v2_preview.png", dpi=300, format="png")
    plt.close(fig)
    print("\nSaved Clean Figure 2 (1000 DPI TIFF & 300 DPI PNG) successfully!")

if __name__ == "__main__":
    main()
