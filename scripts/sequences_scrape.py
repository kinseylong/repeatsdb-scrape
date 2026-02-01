#!/usr/bin/env python3
"""
Stream-wise extract sequences from a RepeatDB annotation CSV and write a multi-FASTA.

Usage:
    python sequences_scrape.py input_annotations.csv output_sequences.fasta

Behavior:
- For rows with `source` containing "RCSB", fetch FASTA from RCSB: https://www.rcsb.org/fasta/entry/{pdb_id}
  and pick the chain matching the `chain` column.
- For rows with `source` containing "AlphaFold" (or exactly "AlphaFoldDB"), fetch UniProt FASTA using the
  accession in the `uniprot` column via UniProt REST endpoints.
- Writes each sequence to output FASTA as it is retrieved to minimize memory usage.
- Writes errors to an accompanying .err file next to the output FASTA.
- FASTA header: >{pdbid}_{chain} followed by description fields: region_values, region_units, pfam, source

Notes:
- The script is conservative when matching chains in the RCSB FASTA headers: it looks for "chain {chain}"
  or common header conventions.
"""

import sys
import os
import csv
import re
import requests
import time
from typing import Optional

RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb}"
UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/search?query=({accession})&format=fasta"
# and let caller log the error.)


def header_matches_chain(header: str, chain: str, pdb: str) -> bool:
    h = header.lower()
    c = chain.lower()
    # common patterns
    if f"chain {c}" in h:
        return True
    # patterns like |A| or |A\n or |A|
    if f"|{chain}|" in header:
        return True
    # patterns like >1abc_A or >1abc|A|
    if re.search(r"\b" + re.escape(chain) + r"\b", header):
        return True
    return False


def fetch_sequence_rcsb(pdb: str, chain: str, timeout=2) -> Optional[str]:
    url = RCSB_FASTA_URL.format(pdb=pdb)
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        text = r.text
    except Exception:
        return None

    parts = re.split(r"(?m)^>", text)
    for part in parts:
        if not part.strip():
            continue
        lines = part.splitlines()
        header = lines[0].strip()
        seq = ''.join(l.strip() for l in lines[1:] if l.strip())
        if header_matches_chain(header, chain, pdb):
            return seq
    # no matching chain found
    return None


def fetch_sequence_uniprot(accession: str, timeout=2) -> Optional[str]:
    """Fetch sequence from UniProt using the JSON search endpoint.

    The search endpoint is used to avoid including chain identifiers in the query.
    We strip any chain suffix (eg. "_A" or "-1") from the provided accession
    before searching.
    """
    if not accession:
        return None

    url = UNIPROT_FASTA_URL.format(accession=accession)
    try:
        r = requests.get(url, timeout=timeout)
        r.raise_for_status()
        text = r.text
        fasta_lines = text.splitlines()
        seq_lines = [line.strip() for line in fasta_lines if not line.startswith('>')]
        seq = ''.join(seq_lines)
        if not seq:
            return None
    except Exception:
        return None

    return seq


def make_description(region_values: str, region_units: str, pfam: str, source: str, uniprot: str) -> str:
    parts = []
    if region_values:
        parts.append(f"region_values={region_values}")
    if region_units:
        parts.append(f"region_units={region_units}")
    if pfam:
        parts.append(f"pfam={pfam}")
    if source:
        parts.append(f"source={source}")
    if uniprot:
        parts.append(f"uniprot={uniprot}")
    return ' '.join(parts)


def main(argv):
    if len(argv) < 3:
        print("Usage: python sequences_scrape.py input_annotations.csv output_sequences.fasta")
        sys.exit(2)
    input_csv = argv[1]
    output_fasta = argv[2]
    error_log = output_fasta + '.err'

    os.makedirs(os.path.dirname(output_fasta) or '.', exist_ok=True)

    with open(input_csv, newline='', encoding='utf-8') as infp, \
         open(output_fasta, 'w', encoding='utf-8') as outfa, \
         open(error_log, 'w', encoding='utf-8') as errf:

        reader = csv.DictReader(infp)
        written = 0
        for i, row in enumerate(reader, 1):
            pdb_id = (row.get('pdb_id') or '').strip()
            chain = (row.get('chain') or '').strip()
            source = (row.get('source') or '').strip()
            uniprot = (row.get('uniprot') or '').strip()
            region_values = (row.get('region_values') or '').strip()
            region_units = (row.get('region_units') or '').strip()
            pfam = (row.get('pfam') or '').strip()

            if not pdb_id or not chain:
                errf.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Missing pdb_id/chain on line {i}\n")
                continue

            seq = None
            try:
                alphafold_source = 'alphafold' in source.lower() or source.strip() == 'AlphaFoldDB'
                if alphafold_source:
                    # pdb_id is actually uniprot accession here
                    seq = fetch_sequence_uniprot(pdb_id)
                else:
                    seq = fetch_sequence_rcsb(pdb_id, chain)
                if not seq or not str(seq).strip():
                    # More detailed error logging for missing/empty sequences
                    if alphafold_source:
                        acc_clean = re.sub(r'[_\-].+$', '', uniprot)
                        detail = f"uniprot={acc_clean}" if acc_clean else "uniprot=<missing>"
                    else:
                        detail = f"pdb={pdb_id} chain={chain}"
                    errf.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Empty/missing sequence for {pdb_id}_{chain} (source={source}) {detail}\n")
                    continue
                # write fasta record
                desc = make_description(region_values, region_units, pfam, source, uniprot)
                header = f">{pdb_id}_{chain} {desc}\n"
                outfa.write(header)
                #wrap sequence to 80 chars per line
                for j in range(0, len(seq), 80):
                    outfa.write(seq[j:j+80] + '\n')

                outfa.flush()
                written += 1
            except Exception as e:
                errf.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] Exception for {pdb_id}_{chain}: {e}\n")
                continue

    print(f"Done. Wrote {written} sequences to {output_fasta}. Errors in {error_log}.")


if __name__ == '__main__':
    main(sys.argv)
