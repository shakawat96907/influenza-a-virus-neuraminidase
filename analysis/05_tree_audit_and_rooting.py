#!/usr/bin/env python3
import re
import pandas as pd
from Bio import Phylo

def parse_iqtree_report(iqtree_file):
    txt = open(iqtree_file).read()
    model = re.search(r"Best-fit model according to BIC:\s*(.+)", txt)
    logl = re.search(r"Log-likelihood of the tree:\s*([-\d.]+)", txt)
    tlen = re.search(r"Total tree length \(sum of branch lengths\):\s*([\d.]+)", txt)
    return {
        "model": model.group(1).strip() if model else "N/A",
        "logl": logl.group(1).strip() if logl else "N/A",
        "tree_length": tlen.group(1).strip() if tlen else "N/A"
    }

def evaluate_support(tree_file):
    tree = Phylo.read(tree_file, "newick")
    internal = [c for c in tree.get_nonterminals() if c.name and "/" in c.name]
    strong_ufboot, strong_both = 0, 0
    for c in internal:
        parts = c.name.split("/")
        try:
            alrt, ufboot = float(parts[0]), float(parts[1])
            if ufboot >= 95: strong_ufboot += 1
            if alrt >= 80 and ufboot >= 95: strong_both += 1
        except (ValueError, IndexError):
            pass
    return tree, len(tree.get_terminals()), len(internal), strong_ufboot, strong_both

def main():
    nuc_rep = parse_iqtree_report("results/trees/tree_nuc_cap20.iqtree")
    prot_rep = parse_iqtree_report("results/trees/tree_prot_cap20.iqtree")

    nuc_tree, n_taxa, n_int, n_uf95, n_both = evaluate_support("results/trees/tree_nuc_cap20.treefile")
    prot_tree, p_taxa, p_int, p_uf95, p_both = evaluate_support("results/trees/tree_prot_cap20.treefile")

    print("=== 1. Global Nucleotide ML Tree (New Figure 1) ===")
    print(f"Taxa count       : {n_taxa}")
    print(f"Best-fit model   : {nuc_rep['model']}")
    print(f"Log-likelihood   : {nuc_rep['logl']}")
    print(f"Total tree length: {nuc_rep['tree_length']}")
    print(f"Well-supported nodes (UFBoot>=95%): {n_uf95}/{n_int} ({n_uf95/max(1,n_int)*100:.1f}%)")
    print(f"Strongly supported (SH-aLRT>=80% & UFBoot>=95%): {n_both}/{n_int} ({n_both/max(1,n_int)*100:.1f}%)\n")

    print("=== 2. Global Amino-Acid ML Tree (Supplementary Figure) ===")
    print(f"Unique taxa count: {p_taxa}")
    print(f"Best-fit model   : {prot_rep['model']}")
    print(f"Log-likelihood   : {prot_rep['logl']}")
    print(f"Total tree length: {prot_rep['tree_length']}")
    print(f"Well-supported nodes (UFBoot>=95%): {p_uf95}/{p_int} ({p_uf95/max(1,p_int)*100:.1f}%)")
    print(f"Strongly supported (SH-aLRT>=80% & UFBoot>=95%): {p_both}/{p_int} ({p_both/max(1,p_int)*100:.1f}%)\n")

    meta = pd.read_csv("data/processed/final_cap20_meta.csv")
    col_map = {c.lower(): c for c in meta.columns}
    acc_col = col_map.get("accession", meta.columns[0])
    date_col = col_map.get("collection_date") or col_map.get("date") or col_map.get("release_date")
    year_col = col_map.get("year")
    country_col = col_map.get("country")

    if year_col:
        meta["_yr"] = pd.to_numeric(meta[year_col], errors="coerce")
    elif date_col:
        meta["_yr"] = pd.to_numeric(meta[date_col].astype(str).str[:4], errors="coerce")
    else:
        meta["_yr"] = 2015

    sort_cols = [c for c in [date_col, acc_col] if c]
    earliest = meta[meta["_yr"] == 2015].sort_values(by=sort_cols).iloc[0]
    root_acc = str(earliest[acc_col]).strip()
    print("=== 3. Outgroup Rooting on Earliest 2015 Isolate ===")
    print(f"Root Accession: {root_acc} | Country: {earliest[country_col] if country_col else 'N/A'} | Date: {earliest[date_col] if date_col else earliest['_yr']}")

    nuc_tree.root_with_outgroup(root_acc)
    Phylo.write(nuc_tree, "results/trees/tree_nuc_cap20_rooted.treefile", "newick")
    print("Saved rooted nucleotide tree to: results/trees/tree_nuc_cap20_rooted.treefile")

if __name__ == "__main__":
    main()
