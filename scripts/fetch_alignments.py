#!/usr/bin/env python3
"""
Lightweight alignment fetcher for RepeatsDB.
Reads an annotations CSV and downloads alignment FASTAs via the RepeatsDB API.
No selenium/beautifulsoup dependencies — only uses requests.
"""

import os
import sys
import csv
import ast
import time
import argparse
import requests
import tqdm


def fetch_alignment(pdb_id, chain, source, region_num, region_id,
                    output_dir="repeatsDB_alignments", error_log=None, silent=False):
    """Download alignment FASTA for a specific region."""
    if source == "AlphaFoldDB":
        db_id = "adb"
        url = (f"https://repeatsdb.org/api/public/production/{db_id}/"
               f"{pdb_id}.{chain}/region.{region_num}/sequence_alignment.fasta")
    else:
        db_id = pdb_id[1:3]
        url = (f"https://repeatsdb.org/api/public/production/pdb/{db_id}/"
               f"{pdb_id}.{chain}/region.{region_num}/sequence_alignment.fasta")

    region_dir_name = f"repeatsDB_alignments_{region_id.replace('.', '_')}"
    region_output_dir = os.path.join(output_dir, region_dir_name)
    os.makedirs(region_output_dir, exist_ok=True)

    output_path = os.path.join(
        region_output_dir, f"{pdb_id}_{chain}_{region_id}_{region_num}.fasta"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        content = r.content
        with open(output_path, "wb") as f:
            f.write(content)
        return content.decode("utf-8", errors="ignore")
    except Exception as e:
        if not silent:
            error_msg = (f"Error fetching {pdb_id} {chain} region {region_id} "
                         f"(region_num={region_num}) from {url}: {e}\n")
            if error_log:
                error_log.write(error_msg)
                error_log.flush()
            else:
                print(error_msg.strip())
        return None


def parse_unit_count(fasta_content):
    """Extract the maximum unit number from FASTA headers like '>unit.5.fasta'."""
    if not fasta_content:
        return 0
    unit_numbers = []
    for line in fasta_content.split("\n"):
        if line.startswith(">unit."):
            parts = line.split(".")
            if len(parts) >= 2:
                try:
                    unit_numbers.append(int(parts[1]))
                except ValueError:
                    continue
    return max(unit_numbers) if unit_numbers else 0


def load_annotations(csv_path):
    """Read annotations CSV and return list of dicts with pdb_id, chain, source, region_values."""
    alignment_data = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pdb = (row.get("pdb_id") or "").strip()
            chain = (row.get("chain") or "").strip()
            source = (row.get("source") or "").strip()
            region_values_str = (row.get("region_values") or "[]").strip()
            try:
                region_values = ast.literal_eval(region_values_str)
            except Exception:
                region_values = []
            if pdb and chain:
                alignment_data.append({
                    "pdb_id": pdb,
                    "chain": chain,
                    "source": source,
                    "region_values": region_values,
                })
    return alignment_data


def main():
    parser = argparse.ArgumentParser(
        description="Download RepeatsDB alignment FASTAs from an annotations CSV"
    )
    parser.add_argument(
        "input_csv",
        help="Path to the annotations CSV file",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to save alignment files (region subdirs created automatically)",
    )
    args = parser.parse_args()

    input_csv = args.input_csv
    output_dir = args.output_dir

    if not os.path.isfile(input_csv):
        print(f"Annotations CSV not found: {input_csv}")
        sys.exit(2)

    os.makedirs(output_dir, exist_ok=True)

    alignment_data = load_annotations(input_csv)
    print(f"Read {len(alignment_data)} annotations from {input_csv}")

    csv_basename = os.path.splitext(os.path.basename(input_csv))[0]
    error_log_path = os.path.join(output_dir, f"{csv_basename}_alignment_errors.log")

    total_downloads = 0
    with open(error_log_path, "a", encoding="utf-8") as error_log:
        error_log.write(f"Started alignments download: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        error_log.write(f"Input CSV: {input_csv}\n")
        error_log.write("=" * 60 + "\n\n")

        for data in tqdm.tqdm(alignment_data, desc=csv_basename):
            pdb_id = data["pdb_id"]
            chain = data["chain"]
            source = data["source"]
            region_values = data.get("region_values", [])

            region_num = 0
            for region_id in region_values:
                if len(region_id) == 1 and region_id.isdigit():
                    continue

                original_region_num = region_num
                fasta_content = None
                max_retries = 30

                for attempt in range(max_retries):
                    current_region_num = original_region_num + attempt
                    fasta_content = fetch_alignment(
                        pdb_id, chain, source,
                        region_num=current_region_num,
                        region_id=region_id,
                        output_dir=output_dir,
                        error_log=error_log,
                        silent=True,
                    )
                    if fasta_content:
                        total_downloads += 1
                        region_num = parse_unit_count(fasta_content)
                        break

                if not fasta_content:
                    error_msg = (
                        f"Error fetching {pdb_id} {chain} region {region_id}: "
                        f"Failed after {max_retries} attempts "
                        f"(region_num {original_region_num} to "
                        f"{original_region_num + max_retries - 1})\n"
                    )
                    error_log.write(error_msg)
                    error_log.flush()
                    region_num = original_region_num

    print(f"Downloaded {total_downloads} alignment files.")
    print(f"Errors logged to: {error_log_path}")


if __name__ == "__main__":
    main()
