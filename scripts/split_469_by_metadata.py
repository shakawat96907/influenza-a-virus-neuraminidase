#!/usr/bin/env python3

import re
import argparse
from Bio import SeqIO

def clean_accession(acc):
    acc = acc.strip()
    acc = acc.replace(">", "")
    acc = acc.split()[0]
    acc = acc.split(".")[0]
    return acc

def extract_accessions_from_line(line):
    """
    Detect accession-like IDs such as BDP39147, QGZ07312, etc.
    """
    pattern = r"\b[A-Z]{2,4}[0-9]{4,}(?:\.[0-9]+)?\b"
    accs = re.findall(pattern, line)
    return [clean_accession(a) for a in accs]

def extract_year_from_line(line):
    """
    Detect year between 2015 and 2025.
    """
    years = re.findall(r"\b(2015|2016|2017|2018|2019|2020|2021|2022|2023|2024|2025)\b", line)
    if years:
        return int(years[0])
    return None

def build_year_map(metadata_file):
    year_map = {}

    with open(metadata_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            year = extract_year_from_line(line)
            accs = extract_accessions_from_line(line)

            if year is not None and accs:
                for acc in accs:
                    year_map[acc] = year

    return year_map

def main():
    parser = argparse.ArgumentParser(description="Split 469 aligned H1N1 NA FASTA using metadata year information.")
    parser.add_argument("--aligned", required=True, help="Aligned FASTA file with accession-only headers")
    parser.add_argument("--metadata", required=True, help="Metadata TXT/CSV file containing accession and year")
    parser.add_argument("--prefix", default="ALN469", help="Output prefix")
    args = parser.parse_args()

    year_map = build_year_map(args.metadata)

    print(f"[INFO] Year mapping found for {len(year_map)} accessions from metadata.")

    out_files = {
        "2015_2017": open(f"{args.prefix}_2015_2017.fasta", "w"),
        "2018_2020": open(f"{args.prefix}_2018_2020.fasta", "w"),
        "2021_2025": open(f"{args.prefix}_2021_2025.fasta", "w"),
        "unassigned": open(f"{args.prefix}_unassigned.fasta", "w")
    }

    counts = {
        "2015_2017": 0,
        "2018_2020": 0,
        "2021_2025": 0,
        "unassigned": 0
    }

    unmapped_examples = []

    for rec in SeqIO.parse(args.aligned, "fasta"):
        acc = clean_accession(rec.id)
        year = year_map.get(acc)

        if year is None:
            SeqIO.write(rec, out_files["unassigned"], "fasta")
            counts["unassigned"] += 1
            if len(unmapped_examples) < 15:
                unmapped_examples.append(acc)

        elif 2015 <= year <= 2017:
            SeqIO.write(rec, out_files["2015_2017"], "fasta")
            counts["2015_2017"] += 1

        elif 2018 <= year <= 2020:
            SeqIO.write(rec, out_files["2018_2020"], "fasta")
            counts["2018_2020"] += 1

        elif 2021 <= year <= 2025:
            SeqIO.write(rec, out_files["2021_2025"], "fasta")
            counts["2021_2025"] += 1

        else:
            SeqIO.write(rec, out_files["unassigned"], "fasta")
            counts["unassigned"] += 1
            if len(unmapped_examples) < 15:
                unmapped_examples.append(acc)

    for f in out_files.values():
        f.close()

    print("\n[SPLIT SUMMARY]")
    total = sum(counts.values())
    for k, v in counts.items():
        print(f"{k}: {v} sequences")
    print(f"Total: {total} sequences")

    if counts["unassigned"] > 0:
        print("\n[WARNING] Some sequences still unassigned.")
        print("Example unmapped accessions:")
        for acc in unmapped_examples:
            print(acc)
        print(f"\nCheck file: {args.prefix}_unassigned.fasta")

if __name__ == "__main__":
    main()
