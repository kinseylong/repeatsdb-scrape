# RepeatsDB Scraper + Sequence Extraction

Small tools to: 1) collect RepeatsDB protein IDs by their classification, 2) download the corresponding MSAs of the repeat-unit, and 3) fetch the corresponding protein sequences.

Database link: 
https://repeatsdb.org/home

## What’s here
- `scripts/repeatsdb_scrape.py` — compiles all proteins within a repeat class, and (optionally) downloads their per-chain repeat-unit MSAs.
- `scripts/sequences_scrape.py` — read an annotations CSV and fetch full-chain sequences from RCSB (PDB) or UniProt/AlphaFold.
- `run_pipeline.py` — run scraper + sequence fetch across region classes. Supports `--regions`, `--save-sequences/--no-save-sequences`, `--save-alignments/--no-save-alignments`, `--workers`, and path overrides.
- `run_pipeline.py` — run scraper + sequence fetch across region classes. Supports `--regions`, `--save-sequences/--no-save-sequences`, `--save-alignments/--no-save-alignments`, `--alignments-only`, `--workers`, and path overrides.

## Outputs
- Annotation CSV: `result-annotations/repeatsDB_annotations_<region>.csv` — RepeatsDB proteins belonging to a certain repeat class (region). Each row corresponds to one protein. Fields include `pdb_id`, `chain`, `source`, `region_values`, `region_units`, `uniprot`, `pfam`, `status`.
- Repeat-unit MSAs: `result-alignments/repeatsDB_alignments_<region>/*` — MSAs of detected repeat units (one FASTA per `pdb_id`_`chain`). These are the aligned repeat segments from RepeatsDB, not full-chain sequences.
- Full-chain sequences: `sequences/repeatsdb_seqs_<region>.fasta` — multi-FASTA of complete chains fetched from PDB or UniProt depending on the ID.
- Logs: `.log`, `_errors.txt`, or `.err` files written alongside outputs for errors and diagnostics.

## Quick start
1. Create a virtualenv and install dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

2. Ensure a compatible Chrome/Chromium is installed for `scripts/repeatsdb_scrape.py` (Selenium + `webdriver-manager`).
	- ChromeDriver must match the installed Chrome version; mismatches commonly cause startup failures. Headless operation is possible but can be fragile (rendering, pagination and JS timing differ).
	- Recommended: run the scraper on a laptop or workstation with Chrome available. For remote/cluster runs, enable headless flags, provide a display (e.g., Xvfb) or a persistent Chrome profile, increase timeouts, and reduce concurrency to avoid failures or rate limits.
    - Note: This tool has only been tested on a macOS with Chrome app installed.


## Usage examples

```bash
# full pipeline for a few region classes
python run_pipeline.py --regions 3.1 3.2 3.3

# scrape one region and save outputs to explicit paths
python scripts/repeatsdb_scrape.py --region-classes 3.3 --output-dir result-alignments/repeatsDB_alignments_3_3 --output-csv result-annotations/repeatsDB_annotations_3_3.csv

# extract full sequences from an existing annotations CSV
python scripts/sequences_scrape.py result-annotations/repeatsDB_annotations_3_3.csv sequences/repeatsdb_seqs_3_3.fasta

# disable saving full sequences
python run_pipeline.py --no-save-sequences

# disable downloading repeat-unit MSAs
python run_pipeline.py --no-save-alignments

# download alignments only for the default regions (requires per-region annotation CSVs to exist)
python run_pipeline.py --alignments-only

# fetch alignments from a single existing annotations CSV (no browser scraping)
python scripts/repeatsdb_scrape.py --output-csv result-annotations/repeatsDB_annotations_3_3.csv --output-dir result-alignments/repeatsDB_alignments_3_3 --fetch-alignments-only
```


## Acknowledgements

Thanks to the RepeatsDB team and dataset contributors for making annotations available; please cite their work when using derived results.

Citation:

Damiano Clementel, Paula Nazarena Arrías, Soroush Mozaffari, et al., "RepeatsDB in 2025: expanding annotations of Structured Tandem Repeats proteins on AlphaFoldDB."

