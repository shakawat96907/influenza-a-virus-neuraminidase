#!/usr/bin/env python3

import re
import argparse
from Bio import SeqIO

def extract_accession(record_id):
    """
    Extract accession from FASTA ID.
    Works for IDs like:
    BDP39147
    gi|xxx|ref|BDP39147|
    sp|BDP39147|...
    """
    rid = record_id.strip()

    if "|" in rid:
        parts = [p for p in rid.split("|") if p.strip()]
        for p in parts:
            if re.match(r"^[A-Za-z]{2,4}[0-9]{4,}", p):
                return p
        return parts[-1]

    return rid.split()[0]

def find_year(text):
    """
    Find year between 2015 and 2025 from FASTA header/description.
    """
    years = re.findall(r"(20[0-2][0-9])", text)
    years = [int(y) for y in years if 2015 <= int(y) <= 2025]

    if len(years) > 0:
        return years[0]
    return None

def build_year_map_from_raw(raw_fasta):
    year_map = {}
    no_year = []

    for rec in SeqIO.parse(raw_fasta, "fasta"):
        acc = extract_accession(rec.id)
        year = find_year(rec.description)

        if year is not None:
            year_map[acc] = year
        else:
            no_year.append(rec.description)

    print(f"[INFO] Year mapping created for {len(year_map)} raw sequences.")
    print(f"[INFO] Raw sequences without year: {len(no_year)}")

    if len(no_year) > 0:
        print("\n[EXAMPLE RAW HEADERS WITHOUT YEAR]")
        for h in no_year[:5]:
            print(h)

    return year_map

def main():
    parser = argparse.ArgumentParser(
        description="Split aligned 469 H1N1 NA protein FASTA into year groups using raw FASTA headers."
    )

    parser.add_argument("--aligned", required=True, help="Aligned 469 FASTA file")
    parser.add_argument("--raw", required=True, help="Raw FASTA file with year in headers")
    parser.add_argument("--prefix", default="ALN469", help="Output prefix")

    args = parser.parse_args()

    year_map = build_year_map_from_raw(args.raw)

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

    unassigned_examples = []

    for rec in SeqIO.parse(args.aligned, "fasta"):
        acc = extract_accession(rec.id)
        year = year_map.get(acc)

        if year is None:
            SeqIO.write(rec, out_files["unassigned"], "fasta")
            counts["unassigned"] += 1
            if len(unassigned_examples) < 10:
                unassigned_examples.append(rec.description)

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
            if len(unassigned_examples) < 10:
                unassigned_examples.append(rec.description)

    for f in out_files.values():
        f.close()

    print("\n[SPLIT SUMMARY]")
    total = sum(counts.values())
    for k, v in counts.items():
        print(f"{k}: {v} sequences")
    print(f"Total: {total} sequences")

    print("\n[OUTPUT FILES]")
    print(f"{args.prefix}_2015_2017.fasta")
    print(f"{args.prefix}_2018_2020.fasta")
    print(f"{args.prefix}_2021_2025.fasta")
    print(f"{args.prefix}_unassigned.fasta")

    if counts["unassigned"] > 0:
        print("\n[WARNING] Some aligned sequences could not be matched with raw FASTA year information.")
        print("Examples:")
        for h in unassigned_examples:
            print(h)

if __name__ == "__main__":
    main()
