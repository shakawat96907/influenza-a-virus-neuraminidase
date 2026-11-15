#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = pd.read_csv("results/tables/Supplementary Table 9 - (regional_enrichment_all_variants).csv")
res["CI_Low"] = res["Fisher_95CI"].str.split("–").str[0].astype(float)
res["CI_High"] = res["Fisher_95CI"].str.split("–").str[1].astype(float)
res["CMH_CI_Low"] = res["CMH_95CI"].str.split("–").str[0].astype(float)
res["CMH_CI_High"] = res["CMH_95CI"].str.split("–").str[1].astype(float)

plot_df = res[res["Global_Pct"] >= 3.0].sort_values("Fisher_OR", ascending=True).tail(18)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.2, 7.4), gridspec_kw={"width_ratios": [1.1, 1.0], "wspace": 0.26})

y_pos = np.arange(len(plot_df))
log_or = np.log2(plot_df["Fisher_OR"])
err_lo = log_or - np.log2(plot_df["CI_Low"])
err_hi = np.log2(plot_df["CI_High"]) - log_or

cmh_log = np.log2(plot_df["CMH_OR"])
cmh_elo = cmh_log - np.log2(plot_df["CMH_CI_Low"])
cmh_ehi = np.log2(plot_df["CMH_CI_High"]) - cmh_log

ax1.axvline(0, color="#333333", linestyle="--", lw=1.0, zorder=1)
h1 = ax1.errorbar(log_or, y_pos + 0.14, xerr=[err_lo, err_hi], fmt="o", color="#C0392B", ecolor="#C0392B", elinewidth=1.3, capsize=2.5, markersize=5.5, label="Unstratified Fisher OR (95% CI)")
h2 = ax1.errorbar(cmh_log, y_pos - 0.14, xerr=[cmh_elo, cmh_ehi], fmt="s", color="#1F618D", ecolor="#1F618D", elinewidth=1.3, capsize=2.5, markersize=5.0, label="Year-Stratified CMH OR (95% CI)")
ax1.set_yticks(y_pos)
ax1.set_yticklabels(plot_df["Variant"], fontsize=9, fontweight="bold")
ax1.set_xlabel("log₂ Odds Ratio (Non-Asia Enriched  ←  0  →  Asia Enriched)", fontsize=9.5, fontweight="bold")
ax1.set_title("A   Regional Enrichment Effect Sizes (Unstratified vs. Year-Stratified)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
ax1.legend(handles=[h1, h2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8, frameon=False)
ax1.grid(axis="x", linestyle=":", alpha=0.4)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

all_v = res.copy()
sig_mask = all_v["Fisher_FDR_P"] < 0.05
s1 = ax2.scatter(all_v.loc[~sig_mask, "NonAsian_Pct"], all_v.loc[~sig_mask, "Asian_Pct"], s=all_v.loc[~sig_mask, "Independent_Origins"] * 0.45 + 18, color="#95A5A6", alpha=0.65, edgecolors="#555555", lw=0.5, label="Non-significant (FDR ≥ 0.05)")
s2 = ax2.scatter(all_v.loc[sig_mask, "NonAsian_Pct"], all_v.loc[sig_mask, "Asian_Pct"], s=all_v.loc[sig_mask, "Independent_Origins"] * 0.45 + 22, color="#D62728", alpha=0.85, edgecolors="#111111", lw=0.7, label="FDR-significant (FDR < 0.05)")
ax2.plot([0, 60], [0, 60], color="#333333", linestyle="--", lw=1.0)

# Custom non-overlapping callouts for key FDR-significant variants
callouts = {
    "A81V": (34.0, 41.5),
    "R77G": (21.0, 39.5),
    "K469N": (39.0, 18.0),
    "V453M": (39.0, 23.5),
    "S339L": (14.5, 4.5),
    "T381I": (3.5, 11.5),
    "V67I": (9.5, 10.5),
    "S450G": (13.5, 0.5),
    "I288V": (9.5, -1.8)
}
for _, r in all_v[sig_mask].iterrows():
    v = r["Variant"]
    if v in callouts:
        tx, ty = callouts[v]
        ax2.annotate(
            v, xy=(r["NonAsian_Pct"], r["Asian_Pct"]), xytext=(tx, ty),
            fontsize=8, fontweight="bold", color="#1A2530",
            arrowprops=dict(arrowstyle="-", lw=0.65, color="#555555", shrinkA=0, shrinkB=2)
        )

ax2.set_xlabel("Non-Asian Control Prevalence (%)", fontsize=9.5, fontweight="bold")
ax2.set_ylabel("Asian Cohort Prevalence (%)", fontsize=9.5, fontweight="bold")
ax2.set_title("B   Regional Prevalence & Recurrence (Size = Indep. Origins)", loc="left", fontsize=10.5, fontweight="bold", pad=26)
ax2.set_xlim(-3, 60)
ax2.set_ylim(-3.5, 60)
ax2.legend(handles=[s1, s2], loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=2, fontsize=8, frameon=False)
ax2.grid(True, linestyle=":", alpha=0.4)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.subplots_adjust(left=0.07, right=0.97, top=0.88, bottom=0.11)
fig.savefig("results/figures/Figure3_Phylogeny_Aware_Regional_Enrichment_v2_1000dpi.tiff", dpi=1000, format="tiff", pil_kwargs={"compression": "tiff_lzw"})
fig.savefig("results/figures/Figure3_Phylogeny_Aware_Regional_Enrichment_v2_preview.png", dpi=300, format="png")
plt.close(fig)
print("Saved clean Figure 3 v2 successfully!")
