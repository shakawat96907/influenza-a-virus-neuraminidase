#!/usr/bin/env python3
"""
Build and verify the FINAL main dataset (cap20): Asian sequences + year-matched
human non-Asian controls.

Run from the project root:   python scripts/build_final_cap20.py

Inputs (read only)
  data/clean_nuc.fasta, data/clean_prot.fasta          first download, QC-passed
  data/clean_nuc_extra.fasta                           second download, QC-passed (control only)
  raw/new/extra_cds_prot.fasta                         optional, used only to cross-check translations
  results/asia_selected_cap20.csv                      Asian selection (757)
  results/control_human_cap20.csv                      human non-Asian control (705)

Outputs
  data/final_cap20_nuc.fasta, data/final_cap20_prot.fasta   FASTA id = accession without version
  data/final_cap20_meta.csv
  results/final_cap20_missing.csv
  results/final_cap20_year_table.csv
  logs/final_cap20_md5.txt
"""
import re, sys, hashlib
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
D, R, RAW, L = ROOT / "data", ROOT / "results", ROOT / "raw" / "new", ROOT / "logs"
L.mkdir(exist_ok=True)

def parse(path):
    """-> {accession_without_version: (version, SeqRecord)}"""
    out = {}
    for rec in SeqIO.parse(path, "fasta"):
        m = re.match(r"lcl\|([A-Z]+\d+)\.(\d+)_(?:cds|prot)_", rec.id)
        if not m:
            raise SystemExit(f"cannot read accession from header: {rec.id}")
        acc = m.group(1)
        if acc in out:
            raise SystemExit(f"duplicate accession in {path}: {acc}")
        out[acc] = (f"{acc}.{m.group(2)}", rec)
    return out

def translate(rec):
    return str(rec.seq.translate(table=1)).rstrip("*")

nuc1 = parse(D / "clean_nuc.fasta")
prot1 = parse(D / "clean_prot.fasta")
nuc2 = parse(D / "clean_nuc_extra.fasta")
assert not (set(nuc1) & set(nuc2)), "first and extra downloads overlap"
assert set(nuc1) == set(prot1), "first-batch nucleotide and protein accession sets differ"

# translation of the first batch must equal the NCBI protein
mism1 = [a for a, (_, r) in nuc1.items() if translate(r) != str(prot1[a][1].seq).rstrip("*")]
print(f"first batch: {len(nuc1)} records, translation == NCBI protein for {len(nuc1) - len(mism1)}, mismatches {len(mism1)}")

# optional cross-check of the extra batch against the raw NCBI protein file
extra_raw = RAW / "extra_cds_prot.fasta"
mism2, checked2 = [], 0
if extra_raw.exists() and extra_raw.stat().st_size > 0:
    p2 = {}
    for rec in SeqIO.parse(extra_raw, "fasta"):
        m = re.match(r"lcl\|([A-Z]+\d+)\.\d+_", rec.id)
        if m: p2[m.group(1)] = str(rec.seq).rstrip("*")
    for a, (_, r) in nuc2.items():
        if a in p2:
            checked2 += 1
            if translate(r) != p2[a]: mism2.append(a)
    print(f"extra batch: cross-checked {checked2} of {len(nuc2)} against raw NCBI protein, mismatches {len(mism2)}")
else:
    print("extra batch: raw/new/extra_cds_prot.fasta not available, protein derived by translation only")

asia = pd.read_csv(R / "asia_selected_cap20.csv", low_memory=False)
ctl = pd.read_csv(R / "control_human_cap20.csv", low_memory=False)
assert asia.Accession.is_unique and ctl.Accession.is_unique and not (set(asia.Accession) & set(ctl.Accession))

asia["group"] = "Asia"
asia["region"] = "asia"
asia["country"] = asia.Geo_Location.astype(str).str.split(":").str[0].str.strip()
asia["year"] = asia.Collection_Date.astype(str).str.extract(r"(\d{4})")[0].astype(int)
ctl["group"] = "control"
ctl = ctl.rename(columns={"Country": "country", "Year": "year"})
cols = ["Accession", "group", "region", "country", "year", "Collection_Date", "Host", "Isolate", "Length", "Submitters"]
meta = pd.concat([asia[cols], ctl[cols]], ignore_index=True)

allnuc = {**nuc1, **nuc2}
miss = meta[~meta.Accession.isin(allnuc)].copy()
miss["reason"] = "not in QC-passed files (removed by QC or not downloaded)"
keep = meta[meta.Accession.isin(allnuc)].copy()

# build outputs
nuc_out, prot_out, ver = [], [], []
for a in keep.Accession:
    v, r = allnuc[a]
    ver.append(v)
    nuc_out.append(SeqRecord(r.seq, id=a, description=""))
    prot_out.append(SeqRecord(Seq(translate(r)), id=a, description=""))
keep["accession_version"] = ver

# independent QC of every kept sequence
bad = []
for a, rn, rp in zip(keep.Accession, nuc_out, prot_out):
    s = str(rn.seq).upper(); p = str(rp.seq)
    amb = sum(ch not in "ACGT" for ch in s) / len(s)
    if not (s.startswith("ATG") and len(s) % 3 == 0 and "*" not in p and amb <= 0.01 and 464 <= len(p) <= 474):
        bad.append(a)

SeqIO.write(nuc_out, D / "final_cap20_nuc.fasta", "fasta")
SeqIO.write(prot_out, D / "final_cap20_prot.fasta", "fasta")
keep.rename(columns={"Accession": "accession", "Collection_Date": "collection_date",
                     "Host": "host", "Isolate": "isolate", "Length": "length_nt",
                     "Submitters": "submitters"}).to_csv(D / "final_cap20_meta.csv", index=False)
miss.rename(columns={"Accession": "accession"}).to_csv(R / "final_cap20_missing.csv", index=False)

# checks
nn = {r.id for r in SeqIO.parse(D / "final_cap20_nuc.fasta", "fasta")}
pp = {r.id: str(r.seq) for r in SeqIO.parse(D / "final_cap20_prot.fasta", "fasta")}
mm = pd.read_csv(D / "final_cap20_meta.csv")
nuc_chk = {r.id: r for r in SeqIO.parse(D / "final_cap20_nuc.fasta", "fasta")}
checks = {
    "(a) nuc and prot have the same accessions": nn == set(pp),
    "(b) accessions equal those in meta": nn == set(mm.accession),
    "(c) no accession twice": len(nn) == len(mm) == len(keep),
    "(d) protein is exact translation of nucleotide": all(translate(nuc_chk[a]) == pp[a] for a in nn),
    "(e) all hosts are Homo sapiens": (mm.host == "Homo sapiens").all(),
    "(f) independent QC passed for every sequence": not bad,
    "(g) no first-batch translation mismatch": not mism1,
    "(h) extra-batch translations equal raw NCBI protein": (not mism2) if checked2 else None,
}
print()
for k, v in checks.items():
    print(("SKIP " if v is None else "PASS " if v else "FAIL ") + k + (" (raw/new/extra_cds_prot.fasta not available)" if v is None else ""))

# year table
tgt_a = asia.groupby("year").size(); tgt_c = ctl.groupby("year").size()
have = mm.groupby(["year", "group"]).size().unstack(fill_value=0)
yt = pd.DataFrame({"asia_target": tgt_a, "asia_final": have.get("Asia"),
                   "control_target": tgt_c, "control_final": have.get("control")}).fillna(0).astype(int)
yt.loc["total"] = yt.sum()
yt.to_csv(R / "final_cap20_year_table.csv")
print("\n", yt.to_string())
print("\nfinal sequences:", len(mm), "| Asia:", int((mm.group == 'Asia').sum()), "| control:", int((mm.group == 'control').sum()))
print("missing from final:", len(miss), "| Asia:", int((miss.group == 'Asia').sum()), "| control:", int((miss.group == 'control').sum()))
print("missing accessions:", miss.Accession.tolist())

with open(L / "final_cap20_md5.txt", "w") as fh:
    for f in ["data/final_cap20_nuc.fasta", "data/final_cap20_prot.fasta", "data/final_cap20_meta.csv"]:
        fh.write(f"{hashlib.md5((ROOT / f).read_bytes()).hexdigest()}  {f}\n")
print("\n" + (L / "final_cap20_md5.txt").read_text())
if any(v is False for v in checks.values()): sys.exit(1)
