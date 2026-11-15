#!/usr/bin/env python3
import os, subprocess
from collections import Counter
import pandas as pd
from Bio import Phylo

def main():
    os.makedirs("results/mugration", exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)

    meta = pd.read_csv("data/processed/final_clean_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col, grp_col, reg_col, yr_col = col_map["accession"], col_map["group"], col_map["region"], col_map["year"]

    reg_map = {
        "asia": "Asia",
        "northamerica": "North_America",
        "europe": "Europe",
        "africa": "Africa",
        "oceania": "Oceania",
        "southamerica": "South_America"
    }
    meta["Region_State"] = meta.apply(
        lambda r: "Asia" if str(r[grp_col]).strip().lower() == "asia" else reg_map.get(str(r[reg_col]).strip().lower(), "Other"),
        axis=1
    )

    states_csv = "results/mugration/region_states.csv"
    meta[[acc_col, "Region_State"]].rename(columns={acc_col: "name", "Region_State": "region"}).to_csv(states_csv, index=False)

    print("[1/3] Running TreeTime mugration on 1,449-taxon global nucleotide tree...")
    subprocess.run([
        "treetime", "mugration",
        "--tree", "results/trees/tree_nuc_clean_rooted.treefile",
        "--states", states_csv,
        "--attribute", "region",
        "--outdir", "results/mugration"
    ], check=True)

    print("\n[2/3] Counting intercontinental transitions...")
    ann_tree = Phylo.read("results/mugration/annotated_tree.nexus", "nexus")

    def get_state(clade):
        if clade.comment:
            for part in clade.comment.strip("[]&").split(","):
                if part.startswith("region="):
                    return part.split("=")[1].strip('"')
        return None

    transitions = Counter()
    asia_to_nonasia, nonasia_to_asia = 0, 0
    for parent in ann_tree.get_nonterminals():
        p_state = get_state(parent)
        if not p_state: continue
        for child in parent.clades:
            c_state = get_state(child)
            if c_state and p_state != c_state:
                transitions[(p_state, c_state)] += 1
                if p_state == "Asia" and c_state != "Asia": asia_to_nonasia += 1
                elif p_state != "Asia" and c_state == "Asia": nonasia_to_asia += 1

    trans_df = pd.DataFrame([{"From_Region": k, "To_Region": k, "Transition_Count": v} for k, v in transitions.most_common()])
    trans_df.to_csv("results/tables/Supplementary Table 6 - (intercontinental_migration_transitions).csv", index=False)

    print("\n=== T5 Intercontinental Migration Summary ===")
    print(f"Total inferred regional transitions : {sum(transitions.values())}")
    print(f"Exports from Asia -> Non-Asia      : {asia_to_nonasia}")
    print(f"Introductions Non-Asia -> Asia     : {nonasia_to_asia}")
    print("\nTop Migration Pathways:")
    print(trans_df.head(10).to_string(index=False))

    print("\n[3/3] Detecting Asia-enriched monophyletic clusters (size >= 10, >= 80% Asia)...")
    state_lookup = dict(zip(meta[acc_col], meta["Region_State"]))
    yr_lookup = dict(zip(meta[acc_col], pd.to_numeric(meta[yr_col], errors="coerce").astype(int)))

    asia_clades = []
    for idx, clade in enumerate(ann_tree.get_nonterminals()):
        tips = [t.name for t in clade.get_terminals()]
        if 10 <= len(tips) <= 300:
            n_asia = sum(1 for t in tips if state_lookup.get(t) == "Asia")
            pct_asia = (n_asia / len(tips)) * 100
            if pct_asia >= 80.0:
                yrs = [yr_lookup[t] for t in tips if t in yr_lookup]
                asia_clades.append({
                    "Clade_ID": f"Asia_Cluster_{idx}",
                    "Total_Tips": len(tips),
                    "Asian_Tips": n_asia,
                    "Asian_Pct": round(pct_asia, 1),
                    "Year_Span": f"{min(yrs)}-{max(yrs)}" if yrs else "N/A"
                })

    clades_df = pd.DataFrame(asia_clades).sort_values(by="Total_Tips", ascending=False)
    clades_df.to_csv("results/tables/Supplementary Table 7 - (asia_enriched_clades).csv", index=False)
    print(f"Identified {len(clades_df)} Asia-enriched clusters (size 10-300, >=80% Asia):")
    if len(clades_df) > 0:
        print(clades_df.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
