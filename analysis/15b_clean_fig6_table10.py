#!/usr/bin/env python3
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

t10 = pd.read_csv("results/tables/Main_Table_10_Glycosylation_Sequons.csv")
t10["Dominant_Motif"] = t10["Dominant_Motif"].apply(lambda s: eval(s) if str(s).startswith("[") else s)

def reclassify(r):
    f1, f2, f3 = r["Pct_2015_2017"], r["Pct_2018_2020"], r["Pct_2021_2025"]
    if f1 >= 80.0 and f2 >= 80.0 and f3 >= 80.0:
        return "Persistent (≥80% across all periods)"
    elif f1 >= 80.0 and f3 < 80.0:
        return "High-frequency declining (≥80% to <80%)"
    elif f1 < 10.0 and f3 >= 50.0:
        return "Gained (<10% to ≥50%)"
    elif f1 >= 50.0 and f3 < 10.0:
        return "Lost (≥50% to <10%)"
    else:
        return "Low-frequency / Sporadic (<10%)"

t10["Temporal_Classification"] = t10.apply(reclassify, axis=1)
t10.to_csv("results/tables/Main_Table_10_Glycosylation_Sequons.csv", index=False)
print(t10.drop(columns=["Observed_Triplets"]).to_string(index=False))

yr_df = pd.read_csv("results/tables/Supplementary Table 12 - (annual_glycosylation_sequons_2015_2025).csv")
years = sorted(yr_df["Year"].unique())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.0, 6.8), gridspec_kw={"width_ratios": [1.05, 1.15], "wspace": 0.25})
x = range(len(t10))
w = 0.26
b1 = ax1.bar([i - w for i in x], t10["Pct_2015_2017"], width=w, color="#2E86C1", label="2015–2017 (n = 497)")
b2 = ax1.bar(x, t10["Pct_2018_2020"], width=w, color="#F39C12", label="2018–2020 (n = 490)")
b3 = ax1.bar([i + w for i in x], t10["Pct_2021_2025"], width=w, color="#C0392B", label="2021–2025 (n = 462)")

ax1.axhline(80, color="#555555", linestyle="--", lw=0.8)
ax1.axhline(50, color="#7F8C8D", linestyle=":", lw=0.8)
ax1.set_xticks(list(x))
ax1.set_xticklabels([f"N{p}" for p in t10["N1_Sequon_Pos"]], rotation=35, ha="right", fontsize=9, fontweight="bold")
ax1.set_ylabel("Sequon Prevalence (%)", fontsize=10, fontweight="bold")
ax1.set_xlabel("Canonical N1 Sequon Position (N-X-S/T, X ≠ P)", fontsize=10, fontweight="bold")
ax1.set_title("A   Predicted N-Linked Glycosylation Sequon Prevalence Across Periods", loc="left", fontsize=10.5, fontweight="bold", pad=38)
ax1.set_ylim(0, 112)
ax1.legend(handles=[b1, b2, b3], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=3, fontsize=8.2, frameon=False)
ax1.grid(axis="y", linestyle=":", alpha=0.4)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

dyn_styles = [
    (42, "#27AE60", "^", "--", "N42 (Persistent)"),
    (50, "#8E44AD", "d", "-", "N50 (Declining)"),
    (386, "#1F618D", "s", "-", "N386 (Low-freq. / Sporadic)"),
    (70, "#E67E22", "o", ":", "N70 (Low-freq. / Sporadic)")
]
handles_b = []
for pos_d, col_d, mk_d, ls_d, lbl_d in dyn_styles:
    sub_y = yr_df[yr_df["N1_Sequon_Pos"] == pos_d].groupby("Year").apply(
        lambda g: (g["Sequon_Count"].sum() / g["Total_N"].sum()) * 100
    )
    ln, = ax2.plot(sub_y.index, sub_y.values, color=col_d, marker=mk_d, linestyle=ls_d, lw=2.0, markersize=5.5, label=lbl_d)
    handles_b.append(ln)

ax2.axvspan(2020.6, 2022.4, color="#F9E79F", alpha=0.4, label="2021–2022 Bottleneck")
ax2.set_xticks(years)
ax2.set_xlabel("Collection Year (2015–2025)", fontsize=10, fontweight="bold")
ax2.set_ylabel("Global Sequon Prevalence (%)", fontsize=10, fontweight="bold")
ax2.set_title("B   Annual Trajectory of Persistent, Declining, and Sporadic Sequons (2015–2025)", loc="left", fontsize=10.5, fontweight="bold", pad=38)
ax2.set_ylim(-4, 108)
ax2.legend(handles=handles_b, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8.0, frameon=False)
ax2.grid(True, linestyle=":", alpha=0.4)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.subplots_adjust(left=0.06, right=0.97, top=0.82, bottom=0.13)
fig.savefig("results/figures/Figure6_Glycosylation_Sequons_v3_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
fig.savefig("results/figures/Figure6_Glycosylation_Sequons_v3_preview.png", dpi=300, format="png")
plt.close(fig)
print("Saved clean Figure 6 v2 successfully!")
