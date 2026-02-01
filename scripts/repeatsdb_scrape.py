import time, re, pandas as pd
import os
import requests
import tqdm
import argparse
from bs4 import BeautifulSoup
import datetime
import sys
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select, WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# ---------- pagination ----------
def get_current_table_html(driver):
    soup = BeautifulSoup(driver.page_source, "html.parser")
    table = soup.find("table")
    return str(table) if table else ""

def set_page_size(driver, value="100", sleep_s=1, timeout=15):
    sel = WebDriverWait(driver, 20).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, "app-pagination-widget select[aria-label='Page size']"))
    )
    prev_table = get_current_table_html(driver)
    Select(sel).select_by_value(value)
    # Wait for table to update; fall back to a short sleep if JS is slow
    try:
        WebDriverWait(driver, timeout).until(lambda d: get_current_table_html(d) != prev_table)
    except Exception:
        # incremental fallback
        total = 0
        waits = [0.5, 1, 2]
        for w in waits:
            time.sleep(w)
            total += w
            if get_current_table_html(driver) != prev_table:
                break

def read_page_index(driver):
    try:
        idx = driver.find_element(By.CSS_SELECTOR, "app-pagination-widget input[formcontrolname='index']")
        return idx.get_attribute("value") or ""
    except:
        return ""

def click_next_page(driver, sleep_s=1):
    """Click next page; return True if navigation likely happened. Uses incremental backoff waiting for page index change."""
    items = driver.find_elements(By.CSS_SELECTOR, "app-pagination-widget ul.pagination li.page-item")
    if not items:
        return False

    # find active li
    active_i = None
    for i, li in enumerate(items):
        if "active" in (li.get_attribute("class") or "").lower():
            active_i = i
            break
    if active_i is None or active_i + 1 >= len(items):
        return False  # no 'next'

    next_li = items[active_i + 1]
    # skip if disabled
    if "disabled" in (next_li.get_attribute("class") or "").lower():
        return False

    before_idx = read_page_index(driver)

    # click
    try:
        clickable = next_li.find_element(By.CSS_SELECTOR, "a,button,span")
    except:
        clickable = next_li
    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", clickable)
    driver.execute_script("arguments[0].click();", clickable)

    # incremental backoff: check for index change, increasing wait if not yet changed
    waits = [sleep_s, sleep_s * 2, sleep_s * 4, sleep_s * 8]
    for w in waits:
        time.sleep(w)
        after_idx = read_page_index(driver)
        if after_idx != before_idx:
            return True
    # final check using WebDriverWait (short)
    try:
        WebDriverWait(driver, 5).until(lambda d: read_page_index(d) != before_idx)
        return True
    except Exception:
        return False

# ---------- scraping ----------
def scrape_annotations(region_classes="3.3", page_size=100, max_pages=None, sleep_s=1, sleep_limit_per_page=20, output_csv=None, profile_dir=None):
    url = f"https://repeatsdb.org/annotations?updated.by=user,predictor,mapping&limit={page_size}&region.classes={region_classes}"
    # Configure Chrome to load faster: disable images/extensions and use eager pageLoadStrategy
    chrome_opts = Options()
    chrome_opts.add_argument("--no-sandbox")
    chrome_opts.add_argument("--disable-dev-shm-usage")
    chrome_opts.add_argument("--disable-extensions")
    chrome_opts.add_argument("--disable-gpu")
    chrome_opts.add_argument("--window-size=1200,800")
    # optional: headless can be enabled if UI is not needed
    # chrome_opts.add_argument("--headless=new")
    # block images to speed up network/load
    chrome_prefs = {"profile.managed_default_content_settings.images": 2}
    chrome_opts.add_experimental_option("prefs", chrome_prefs)
    # prefer DOM readiness over full resource load
    chrome_opts.set_capability("pageLoadStrategy", "eager")

    # Reuse a persistent Chrome profile if requested to avoid startup overhead
    if profile_dir:
        profile_dir = os.path.expanduser(profile_dir)
        profile_dir = os.path.abspath(profile_dir)
        os.makedirs(profile_dir, exist_ok=True)
        chrome_opts.add_argument(f"--user-data-dir={profile_dir}")

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=chrome_opts)
    
    # Track unique records and minimal data for alignments
    seen_keys = set()
    alignment_data = []  # (pdb_id, chain, source)
    total_records = 0
    
    # Open CSV file for writing
    csv_file = None
    csv_writer = None
    if output_csv:
        csv_file = open(output_csv, 'w', newline='', encoding='utf-8')
        # Write header
        csv_file.write("index,pdb_id,chain,source,region_values,region_units,uniprot,pfam,status\n")
    
    try:
        driver.get(url)
        # Wait for table or at least its tbody rows to appear (faster than sleeping)
        try:
            WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.CSS_SELECTOR, "table tbody tr")))
        except Exception:
            # fallback short sleep if waiting failed
            time.sleep(sleep_s)
        # Capture initial page results (they appear immediately on first load)
        initial_html = get_current_table_html(driver)
        initial_records = parse_table(initial_html)
        for rec in initial_records:
            key = (rec['pdb_id'], rec['chain'])
            if key not in seen_keys:
                seen_keys.add(key)
                total_records += 1
                alignment_data.append({
                    'pdb_id': rec['pdb_id'],
                    'chain': rec['chain'],
                    'source': rec['source']
                })
                if csv_file:
                    region_values_str = str(rec['region_values'])
                    region_units_str = str(rec['region_units'])
                    pfam_str = str(rec['pfam'])
                    csv_file.write(f'"{rec["index"]}","{rec["pdb_id"]}","{rec["chain"]}","{rec["source"]}","{region_values_str}","{region_units_str}","{rec["uniprot"]}","{pfam_str}","{rec["status"]}"\n')

        # Then set page size (may reload table); dedup via seen_keys
        set_page_size(driver, str(page_size), sleep_s=sleep_s)

        page_num = 0
        while max_pages is None or page_num < max_pages:
            total_sleep_time_per_page = 0
            print("Scraping page ", page_num + 1, "...")
            # scrape current page
            html = get_current_table_html(driver)
            page_records = parse_table(html)

            # If no records, wait and try again (up to sleep limit)
            while not page_records:
                    time.sleep(sleep_s)
                    total_sleep_time_per_page += sleep_s
                    new_html = get_current_table_html(driver)
                    page_records = parse_table(new_html)
                    if page_records or total_sleep_time_per_page >= sleep_limit_per_page:
                        print(f"Waited {total_sleep_time_per_page}s for page {page_num + 1}, continuing...")
                        break

            if not page_records:
                print(f"No records found on page {page_num + 1}, ending scrape.")
                break

            # Write records to CSV and track for alignments
            for rec in page_records:
                key = (rec['pdb_id'], rec['chain'])
                if key not in seen_keys:
                    seen_keys.add(key)
                    total_records += 1

                    # Save minimal data for alignments
                    alignment_data.append({
                        'pdb_id': rec['pdb_id'],
                        'chain': rec['chain'],
                        'source': rec['source']
                    })

                    # Write to CSV
                    if csv_file:
                        # Convert lists to string representation
                        region_values_str = str(rec['region_values'])
                        region_units_str = str(rec['region_units'])
                        pfam_str = str(rec['pfam'])

                        csv_file.write(f'"{rec["index"]}","{rec["pdb_id"]}","{rec["chain"]}","{rec["source"]}","{region_values_str}","{region_units_str}","{rec["uniprot"]}","{pfam_str}","{rec["status"]}"\n')

            # If fewer records than page_size, assume this is the last page
            per_page_count = len(page_records)
            if per_page_count < page_size:
                break

            # If we have exactly page_size results, attempt to go to next page
            if per_page_count == page_size:
                # try to go to next page; sleep 3s in click_next_page
                if not click_next_page(driver, sleep_s=sleep_s):
                    break

                page_num += 1
                continue

            # Default: stop if none of the above conditions matched
            break

        return alignment_data, total_records
    finally:
        driver.quit()
        if csv_file:
            csv_file.close()

# ---------- table parsing ----------
NUM_RE = re.compile(r"\d+(?:\.\d+)*(?!\s*units)", re.I)
UNITS_RE = re.compile(r"\b\d+\s*units\b", re.I)

def parse_row(tr):
    tds = tr.find_all("td", recursive=False)
    if len(tds) < 8:
        return None

    index = tds[0].get_text(strip=True)
    img = tds[1].find("img")
    #preview_img = img["src"] if img and img.has_attr("src") else None
    pdb_id = tds[2].get_text(strip=True)
    chain = tds[3].get_text(strip=True)
    source = tds[4].get_text(" ", strip=True)

    # multiple region spans
    region_cell = tds[5]
    region_values, region_units = [], []
    for region_span in region_cell.select(".text-bg-region"):
        region_text = region_span.get_text(" ", strip=True)
        for m in NUM_RE.findall(region_text):
            region_values.append(m)  # keep as string "3.3.1" etc
        for m in UNITS_RE.findall(region_text):
            region_units.append(m)

    ext_cell = tds[6]
    badges = [b.get_text(" ", strip=True) for b in ext_cell.select(".badge")]
    uniprot = None
    pfam = []
    for b in badges:
        if "UniProt" in b:
            uniprot = b.split()[0]
        elif "Pfam" in b:
            pfam.append(b.split()[0])

    status = tds[7].get_text(strip=True)

    return {
        "index": index,
        #"preview_img": preview_img,
        "pdb_id": pdb_id,
        "chain": chain,
        "source": source,
        "region_values": region_values,  # e.g. ["4.4", "3.3.1"]
        "region_units": region_units,    # e.g. ["6 units", "3 units"]
        "uniprot": uniprot,
        "pfam": pfam,
        "status": status,
    }

def parse_table(table_html):
    soup = BeautifulSoup(table_html, "html.parser")
    rows = soup.select("tbody tr") or soup.select("tr")
    out = []
    for tr in rows:
        if tr.select_one("td[rowspan] img[src*='preview']") or tr.select_one("td img[src*='preview']"):
            rec = parse_row(tr)
            if rec: out.append(rec)
    return out


# ---------- download MSAs ----------
def fetch_alignment(pdb_id, chain, source, output_dir="repeatsDB_alignments", error_log=None):
    """
    Example download for the alignment FASTA. Adjust endpoint if needed.
    Directory shard is pdb_id[1:3], e.g. '1a17' -> 'a1' as seen in preview path.
    """
    if source == "AlphaFoldDB":
        db_id = "adb"
        url = f"https://repeatsdb.org/api/public/production/{db_id}/{pdb_id}.{chain}/region.0/sequence_alignment.fasta"
    else:
        db_id = pdb_id[1:3]  # 'a1' for '1a17'
        url = f"https://repeatsdb.org/api/public/production/pdb/{db_id}/{pdb_id}.{chain}/region.0/sequence_alignment.fasta"

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{pdb_id}_{chain}.fasta")
    
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(output_path, "wb") as f:
            f.write(r.content)
    except Exception as e:
        error_msg = f"Error fetching {pdb_id} {chain} from {url}: {e}\n"
        if error_log:
            error_log.write(error_msg)
            error_log.flush()  # Ensure it's written immediately
        else:
            print(error_msg.strip())


# ---------- main ----------
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Scrape RepeatsDB annotations and alignments")
    parser.add_argument(
        "--region-classes",
        type=str,
        default="3.3",
        help="Region classes to filter by (default: 3.3)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory to save alignment files (default: result-alignments/repeatsDB_alignments_YYYYMMDD)"
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default=None,
        help="Path to save CSV output (default: result-annotations/repeatsdb_annotations_YYYYMMDD.csv)"
    )
    parser.add_argument(
        "--chrome-profile-dir",
        type=str,
        default=None,
        help="Path to Chrome user-data-dir to reuse profile and cache (optional)"
    )
    parser.add_argument(
        "--max-pages",
        type=int,
        default=None,
        help="Maximum number of pages to scrape (default: None = scrape all pages)"
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=100,
        help="Number of results per page (default: 100)"
    )
    parser.add_argument(
        "--skip-alignments",
        action="store_true",
        help="Skip downloading alignment files"
    )
    
    args = parser.parse_args()
    
    ### inputs
    page_size = args.page_size
    max_pages = args.max_pages  # Keep as None if not specified
    get_alignments = not args.skip_alignments
    region_classes = args.region_classes
    profile_dir = args.chrome_profile_dir

    job_id = datetime.datetime.now().strftime("%Y%m%d")
    output_dir = args.output_dir or f"./result-alignments/repeatsDB_alignments_{job_id}"
    output_csv = args.output_csv or f"./result-annotations/repeatsDB_annotations_{job_id}.csv"
    output_csv = os.path.expanduser(output_csv)  # Expand ~ and make absolute
    output_dir = os.path.expanduser(output_dir)
    
    # Create directory for CSV if it doesn't exist
    csv_dir = os.path.dirname(output_csv)
    if csv_dir:
        os.makedirs(csv_dir, exist_ok=True)

    # Prepare combined error/log file (append so existing file is preserved)
    error_log_path = output_csv.replace('.csv', '.log')

    # Tee class to write to both console and log file
    class Tee:
        def __init__(self, *writers):
            self.writers = writers
        def write(self, data):
            for w in self.writers:
                try:
                    w.write(data)
                except Exception:
                    pass
            try:
                for w in self.writers:
                    if hasattr(w, 'flush'):
                        w.flush()
            except Exception:
                pass
        def flush(self):
            for w in self.writers:
                try:
                    if hasattr(w, 'flush'):
                        w.flush()
                except Exception:
                    pass

    # Open (or create) the combined errors/log file and redirect stdout/stderr
    log_file = open(error_log_path, 'a', encoding='utf-8')
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = Tee(original_stdout, log_file)
    sys.stderr = Tee(original_stderr, log_file)

    print("Starting scrape...")

    try:
        # Write a header to indicate a new run was appended
        print("\n" + "="*60)
        print(f"New run: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Region classes: {region_classes}")

        alignment_data, total_records = scrape_annotations(
            region_classes=region_classes, 
            page_size=page_size, 
            max_pages=max_pages, 
            sleep_s=5,
            output_csv=output_csv,
            profile_dir=profile_dir
        )
    finally:
        # allow later parts to reopen the file as needed; keep log_file open until end
        pass

    print(f"Scraped {total_records} unique proteins.")

    if get_alignments:
        print(f"Downloading alignments...")
        # append to existing error log instead of overwriting so all prints are combined
        with open(error_log_path, 'a', encoding='utf-8') as error_log:
            error_log.write(f"Started alignments download: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            error_log.write("="*60 + "\n\n")

            for data in tqdm.tqdm(alignment_data, total=len(alignment_data)):
                #check if alignment is already downloaded
                if not os.path.isfile(os.path.join(output_dir, f"{data['pdb_id']}_{data['chain']}.fasta")):
                    fetch_alignment(data["pdb_id"], data["chain"], data["source"], output_dir=output_dir, error_log=error_log)

        print(f"Downloaded alignments for {len(alignment_data)} proteins.")
        print(f"Errors logged to: {error_log_path}")

    print(f"Job complete.")

    # restore stdout/stderr and close the log file
    try:
        sys.stdout = original_stdout
        sys.stderr = original_stderr
    except Exception:
        pass
    try:
        log_file.close()
    except Exception:
        pass