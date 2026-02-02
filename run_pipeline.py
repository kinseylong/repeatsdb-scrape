#!/usr/bin/env python3
"""
Pipeline runner: run repeatdb scraping and sequence extraction for multiple region classes.

Usage:
    python run_pipeline.py [--regions 3.1 3.2 ...] [--no-sequences] [--no-alignments]

Defaults:
  regions: ["3.1","3.2","3.3","3.4","4.1","4.2","4.3","4.4","4.5","4.6","4.7","4.8","5.1","5.2","5.3","5.4","5.5"]
  save sequences: True -> directory `sequences/`
  save alignments: True -> directory `result-alignments/`

Assumes helper scripts (`repeatsdb_scrape.py`, `sequences_scrape.py`) are in `scripts/` by default; falls back
to project root if not present.
"""

import argparse
import os
import subprocess
import sys
import concurrent.futures
import multiprocessing
from datetime import datetime
from typing import List

DEFAULT_REGIONS = [
    "3.1", "3.2", "3.3", "3.4",
    "4.1", "4.2", "4.3", "4.4", "4.5", "4.6", "4.7", "4.8",
    "5.1", "5.2", "5.3", "5.4", "5.5"
]


# helper to normalize region string for filenames (3.1 -> 3_1)
def region_fname(region: str) -> str:
    return region.replace('.', '_')


def run_region(region: str,
               annotations_dir: str,
               alignments_dir: str,
               sequences_dir: str,
               scripts_dir: str,
               save_alignments: bool,
               save_sequences: bool,
               python_executable: str,
               alignments_only: bool = False) -> int:
    """Run repeatsdb_scrape for `region`, then run sequences_scrape on its CSV.
    Returns 0 on success for both steps, non-zero if repeatsdb_scrape fails.
    """
    region_id = region_fname(region)

    # paths per-region
    annotation_csv = os.path.join(annotations_dir, f"repeatsDB_annotations_{region_id}.csv")
    alignment_outdir = os.path.join(alignments_dir, f"repeatsDB_alignments_{region_id}")
    sequences_out = os.path.join(sequences_dir, f"repeatsDB_seqs_{region_id}.fasta")

    os.makedirs(os.path.dirname(annotation_csv) or '.', exist_ok=True)
    os.makedirs(alignment_outdir, exist_ok=True)

    python = python_executable or 'python'

    # build repeatsdb_scrape command
    if alignments_only:
        scrape_cmd = [python, os.path.join(scripts_dir, 'repeatsdb_scrape.py'),
                      '--output-dir', alignment_outdir,
                      '--output-csv', annotation_csv,
                      '--fetch-alignments-only']
    else:
        scrape_cmd = [python, os.path.join(scripts_dir, 'repeatsdb_scrape.py'),
                      '--region-classes', region,
                      '--output-dir', alignment_outdir,
                      '--output-csv', annotation_csv]

        if not save_alignments:
            scrape_cmd.append('--skip-alignments')

    print(f"Running repeatsdb_scrape for region {region} -> {annotation_csv}")
    r = subprocess.run(scrape_cmd)
    if r.returncode != 0:
        print(f"repeatsdb_scrape failed for region {region} with exit {r.returncode}")
        return r.returncode

    # If alignments-only mode was requested, skip sequence extraction
    if alignments_only:
        print(f"Alignments-only mode: skipping sequence extraction for region {region}.")
        return 0

    if save_sequences:
        # run sequences_scrape on the produced CSV
        seq_cmd = [python, os.path.join(scripts_dir, 'sequences_scrape.py'), annotation_csv, sequences_out]
        print(f"Running sequences_scrape for region {region} -> {sequences_out}")
        s = subprocess.run(seq_cmd)
        if s.returncode != 0:
            print(f"sequences_scrape failed for region {region} with exit {s.returncode}")
            # do not treat as fatal for other regions
    else:
        print(f"Skipping sequences_scrape for region {region} (--no-sequences)")

    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="Run RepeatDB pipeline: scrape + extract sequences")
    parser.add_argument("--regions", nargs="*", default=DEFAULT_REGIONS,
                        help="List of region classes to process (default: common list)")
    parser.add_argument("--no-sequences", action="store_true", help="Skip extracting sequences")
    parser.add_argument("--no-alignments", action="store_true", help="Skip downloading alignments")
    parser.add_argument("--alignments-only", action="store_true",
                        help="Do not run the web scraper; read existing annotation CSVs and download alignments only")
    parser.add_argument("--sequences-dir", type=str, default="sequences", help="Directory for FASTA outputs")
    parser.add_argument("--alignments-dir", type=str, default="result-alignments", help="Directory for alignment outputs")
    parser.add_argument("--annotations-dir", type=str, default="result-annotations", help="Directory to write/read annotation CSVs")
    parser.add_argument("--scripts-dir", type=str, default="scripts", help="Directory containing helper scripts")
    parser.add_argument("--workers", type=int, default=1, help="Number of regions to process in parallel (default: 1)")

    args = parser.parse_args(argv)

    regions = args.regions
    save_sequences = not args.no_sequences
    save_alignments = not args.no_alignments
    alignments_only = args.alignments_only
    sequences_dir = args.sequences_dir
    alignments_dir = args.alignments_dir
    annotations_dir = args.annotations_dir
    scripts_dir = args.scripts_dir

    os.makedirs(sequences_dir, exist_ok=True)
    os.makedirs(alignments_dir, exist_ok=True)
    os.makedirs(annotations_dir, exist_ok=True)

    # concurrency setup
    workers = max(1, min(args.workers, len(regions))) if hasattr(args, 'workers') else 1

    print(f"Processing {len(regions)} regions with {workers} worker(s)")

    python_exec = sys.executable or 'python'

    if workers == 1:
        for region in regions:
            run_region(region, annotations_dir, alignments_dir, sequences_dir, scripts_dir,
                       save_alignments, save_sequences, python_exec, alignments_only=alignments_only)
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as exe:
            futures = {
                exe.submit(run_region, region, annotations_dir, alignments_dir, sequences_dir, scripts_dir,
                           save_alignments, save_sequences, python_exec, alignments_only): region
                for region in regions
            }
            for fut in concurrent.futures.as_completed(futures):
                region = futures[fut]
                try:
                    code = fut.result()
                    if code != 0:
                        print(f"Region {region} finished with non-zero code {code}")
                except Exception as e:
                    print(f"Region {region} raised exception: {e}")


if __name__ == "__main__":
    main()