# RepeatDB Scrape + Sequence Extraction

This repository contains small scripts to scrape RepeatsDB annotations and download associated alignments and protein sequences.

Contents
- `repeatdb_scrape.py` — Scrapes RepeatsDB web UI using Selenium to collect annotation rows and (optionally) download per-chain alignments. Writes per-run CSV and alignment FASTA files. Defaults to writing outputs under `annotations/` and `alignments/` paths with a date-based job id.
- `sequences_scrape.py` — Reads an annotations CSV (produced by `repeatdb_scrape.py`) and fetches protein sequences for each `pdb_id,chain` pair using RCSB (for PDB entries) or UniProt (for AlphaFold/UniProt entries). Streams FASTA output and logs errors to a `.err` file.
- `run_all_regions.py` — Simple helper to run `repeatdb_scrape.py` for a fixed list of region classes.
- `run_bulk_sequences_scrape.py` — Runs `sequences_scrape.py` on all `repeatdb_annotations_*.csv` files inside `result-annotations/`, writing FASTA outputs to `sequences/`.
- `run_pipeline.py` — Pipeline runner that coordinates scraping and sequence extraction for multiple region classes. Supports CLI args to control which steps to run and where outputs go.

Quick notes on outputs
- **Annotation CSVs**: `result-annotations/repeatsDB_annotations_<region>.csv` — CSVs produced by `repeatsdb_scrape.py` that contain RepeatsDB annotation rows (one row per detected repeat instance). Typical fields include `pdb_id`, `chain`, `source`, `region_values`, `region_units`, `uniprot`, `pfam`, and `status`.
- **Alignment FASTAs (MSAs)**: `result-alignments/repeatsDB_alignments_<region>/*` — These are the multiple-sequence alignments (MSAs) of the detected repeat units downloaded from RepeatsDB. Files are written per `pdb_id`_`chain` (one `.fasta` per chain) and represent the aligned repeat segments, not the full protein sequence.
- **Sequence FASTA (full protein sequences)**: `sequences/repeatsdb_seqs_<region>.fasta` — Multi-sequence FASTA files written by `sequences_scrape.py`. These contain the full protein sequences (complete chains) retrieved from PDB (via RCSB) or from UniProt/AlphaFold entries, and are distinct from the repeat MSAs above.
- **Error logs**: next to outputs with suffixes like `_errors.txt`, `.log`, or `.err` for sequence fetch errors and alignment download problems.

Why this needs a laptop / workstation (limitations)
- `repeatdb_scrape.py` uses Selenium + ChromeDriver and interacts with the RepeatsDB web UI. This typically requires:
  - A real Chrome/Chromium binary installed on the machine, and
  - A display (X11/Quartz) or a working headless Chrome configuration.

Because the script drives the real browser and relies on page rendering and JS-driven pagination, it is fragile on HPC clusters that do not provide a display or do not allow launching browser processes. The `webdriver-manager` package is used to fetch a matching ChromeDriver automatically, but you still need Chrome installed and available in PATH.

If you plan to run scraping on a remote server or cluster, consider one of these options:
- Run the scraper locally on your laptop (recommended) where a GUI and Chrome are available.
- Configure Chrome to run headless (modify `repeatdb_scrape.py` to add ChromeOptions with `--headless=new` or `--headless` and appropriate flags). This has not been tested: headless mode is not guaranteed to behave identically to a real desktop browser.

Setup (recommended)
1. Create and activate a virtual environment (venv):

```bash
python -m venv .venv
source .venv/bin/activate
```

2. Install dependencies (the project may include a `requirements.txt`):

```bash
pip install -r requirements.txt
```

3. Ensure Chrome/Chromium is installed and up-to-date. On macOS you can install Chrome from Google or via Homebrew:

```bash
brew install --cask google-chrome
```

4. (Optional) If you need headless operation, edit `repeatdb_scrape.py` to add ChromeOptions. Example:

```python
from selenium.webdriver.chrome.options import Options
opts = Options()
opts.add_argument('--headless=new')
driver = webdriver.Chrome(options=opts, service=Service(ChromeDriverManager().install()))
```

Usage examples
- Scrape annotations + alignments for region class `3.3` (default job id uses YYYYMMDD):

```bash
python repeatdb_scrape.py --region-classes 3.3 --output-dir result-alignments/repeatDB_alignments_20260131 --output-csv result-annotations/repeatdb_annotations_3_3.csv
```

- Extract sequences from a CSV and write FASTA as you go:

```bash
python sequences_scrape.py result-annotations/repeatdb_annotations_3_3.csv sequences/repeatsdb_seqs_3_3.fasta
```

- Run the full pipeline (scrape then sequences) for the default region set:

```bash
python run_pipeline.py
```

- Run only sequence extraction from existing CSVs:

```bash
python run_pipeline.py --no-alignments
```

Notes & troubleshooting
- If you see an urllib3/OpenSSL/LibreSSL warning, it relates to the Python environment's SSL backend. It is usually a warning and doesn't block HTTP requests, but consider using a system Python or conda environment with OpenSSL if you see TLS errors.
- If ChromeDriver fails to start, check that Chrome is installed and compatible with the driver the `webdriver-manager` installs.