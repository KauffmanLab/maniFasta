#!/usr/bin/env python3

"""
mod_II_fetchProteinAccessions_exact.py
───────────────────────────────
Fetch protein FASTA records from NCBI, rewrite each FASTA first token to the
exact maniFasta protein_id, and write an explicit provenance report.

This helper intentionally does no alias matching and no identifier guessing.
For each input row, the requested fetch_id must return exactly one FASTA record.
The output FASTA first token is always the row's protein_id.

To keep exact one-row → one-record provenance while avoiding slow serial network
latency, records are fetched concurrently with a shared request pacer. Each
worker fetches exactly one input row by its own fetch_id; outputs are written in
input order after all exact fetches complete.

Input TSV columns:
  Required fallback:
    accession

  Optional explicit provenance columns:
    protein_id        ID to write as FASTA first token and metadata join key
    fetch_id          ID requested from NCBI; defaults to accession
    fetch_type        e.g. GI, accession; defaults to accession
    source_accession  original source accession/name; defaults to accession
    source_label      source label for report/logging

Output:
  --output FASTA with headers rewritten as:
      >protein_id fetch_id=... fetch_type=... source_accession=... returned_fasta_id=... <returned description>

  --report TSV containing one row per requested record.
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import csv
import sys
import threading
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOLNAME = "maniFasta_fetchProteinExact"


def info(msg: str) -> None:
    print(f"[INFO] {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"[WARN] {msg}", file=sys.stderr)


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def clean_token(value: str) -> str:
    """Return a FASTA-header-safe key/value token value with no whitespace."""
    value = (value or "NA").strip()
    if not value:
        value = "NA"
    return "_".join(value.replace("\r", " ").replace("\n", " ").split())


def clean_desc(value: str) -> str:
    return " ".join((value or "").replace("\r", " ").replace("\n", " ").split())


def read_input(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader((line for line in fh if line.strip() and not line.startswith("#")), delimiter="\t")
        if reader.fieldnames is None:
            die(f"{path} has no header row")
        lower_to_name = {h.strip().lower(): h for h in reader.fieldnames}
        if "accession" not in lower_to_name and "fetch_id" not in lower_to_name:
            die(f"{path} must contain at least an accession or fetch_id column. Found: {reader.fieldnames}")
        rows = []
        for n, row in enumerate(reader, start=2):
            r = {k: (v or "").strip().replace("\r", "") for k, v in row.items()}

            def get(*names: str) -> str:
                for name in names:
                    key = lower_to_name.get(name.lower())
                    if key is not None and r.get(key):
                        return r[key].strip()
                return ""

            accession = get("accession")
            fetch_id = get("fetch_id") or accession
            protein_id = get("protein_id") or fetch_id or accession
            fetch_type = get("fetch_type") or ("accession" if accession else "NA")
            source_accession = get("source_accession") or accession or fetch_id
            source_label = get("source_label") or "NA"
            if not fetch_id:
                die(f"row {n}: empty fetch_id/accession")
            if not protein_id:
                die(f"row {n}: empty protein_id")
            rows.append({
                "source_label": source_label,
                "protein_id": protein_id,
                "fetch_id": fetch_id,
                "fetch_type": fetch_type,
                "source_accession": source_accession,
                "accession": accession or fetch_id,
                "input_row": str(n),
            })
    return rows


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header = None
    seq_chunks: List[str] = []
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_chunks)))
            header = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_chunks)))
    return records


class RequestPacer:
    """Thread-safe pacer that spaces request starts by a fixed interval."""

    def __init__(self, interval_sec: float):
        self.interval_sec = max(0.0, float(interval_sec))
        self.lock = threading.Lock()
        self.next_allowed = 0.0

    def wait(self) -> None:
        if self.interval_sec <= 0:
            return
        with self.lock:
            now = time.monotonic()
            wait_sec = max(0.0, self.next_allowed - now)
            self.next_allowed = max(now, self.next_allowed) + self.interval_sec
        if wait_sec > 0:
            time.sleep(wait_sec)


BATCH_SIZE = 200  # NCBI efetch POST limit for protein db


def fetch_batch(fetch_ids: list[str], email: str | None, api_key: str | None,
                retries: int, pacer: RequestPacer) -> str:
    """Fetch up to BATCH_SIZE protein records in a single POST. Returns raw FASTA text."""
    params = {
        "db":      "protein",
        "id":      ",".join(fetch_ids),
        "rettype": "fasta",
        "retmode": "text",
        "tool":    TOOLNAME,
    }
    if email:   params["email"]   = email
    if api_key: params["api_key"] = api_key
    data = urllib.parse.urlencode(params).encode("utf-8")
    url  = f"{EUTILS}/efetch.fcgi"

    last_err = None
    for attempt in range(1, retries + 1):
        pacer.wait()
        try:
            req = urllib.request.Request(url, data=data, method="POST")
            with urllib.request.urlopen(req, timeout=120) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as exc:
            last_err = exc
            warn(f"batch attempt {attempt}/{retries} failed ({len(fetch_ids)} ids): {exc}")
            if attempt < retries:
                time.sleep(2 ** attempt)
    raise RuntimeError(f"batch fetch failed after {retries} attempts: {last_err}")


def wrap(seq: str, width: int = 60) -> Iterable[str]:
    for i in range(0, len(seq), width):
        yield seq[i:i+width]


def process_batch(batch: list[tuple[int, Dict[str, str]]], *,
                  email: str | None, api_key: str | None,
                  retries: int, pacer: RequestPacer) -> list[Dict[str, object]]:
    """Fetch a batch of rows in a single NCBI POST and map each returned record
    back to its input row by fetch_id (with version-stripped fallback).

    Returns result dicts in the same order as the input batch.
    """
    fetch_ids = [row["fetch_id"] for _, row in batch]
    text = fetch_batch(fetch_ids, email, api_key, retries, pacer)
    records = parse_fasta(text)

    # Index returned records by their first FASTA token for O(1) lookup.
    returned: Dict[str, tuple[str, str]] = {}
    for hdr, seq in records:
        rid = hdr.split()[0] if hdr.split() else ""
        if rid:
            returned[rid] = (hdr, seq)

    results = []
    for idx, row in batch:
        protein_id     = row["protein_id"]
        fetch_id       = row["fetch_id"]
        fetch_type     = row["fetch_type"]
        source_acc     = row["source_accession"]

        # Match: exact fetch_id first, then version-stripped (e.g. NP_001234.2 → NP_001234).
        matched_hdr, matched_seq = returned.get(fetch_id, (None, None))
        if matched_hdr is None:
            bare = fetch_id.split(".")[0]
            for rid, (hdr, seq) in returned.items():
                if rid.split(".")[0] == bare:
                    matched_hdr, matched_seq = hdr, seq
                    break

        if matched_hdr is None or not matched_seq:
            raise RuntimeError(
                f"fetch_id={fetch_id!r} not found in batch response "
                f"({len(fetch_ids)} requested, {len(records)} returned)"
            )

        parts        = matched_hdr.split()
        returned_id  = parts[0] if parts else "NA"
        returned_desc = " ".join(parts[1:]) if len(parts) > 1 else ""

        out_header = (
            f">{clean_token(protein_id)} "
            f"fetch_id={clean_token(fetch_id)} "
            f"fetch_type={clean_token(fetch_type)} "
            f"source_accession={clean_token(source_acc)} "
            f"returned_fasta_id={clean_token(returned_id)}"
        )
        desc = clean_desc(returned_desc)
        if desc:
            out_header += f" {desc}"

        results.append({
            "idx":    idx,
            "header": out_header,
            "seq":    matched_seq,
            "report_row": {
                "source_label":         row["source_label"],
                "protein_id":           protein_id,
                "fetch_id":             fetch_id,
                "fetch_type":           fetch_type,
                "source_accession":     source_acc,
                "returned_fasta_id":    returned_id,
                "returned_description": returned_desc,
                "sequence_length":      len(matched_seq),
                "status":               "OK",
                "note":                 "",
            },
        })

    return results


def main() -> None:
    ap = argparse.ArgumentParser(description="Fetch NCBI protein records with exact protein_id rewriting and provenance.")
    ap.add_argument("-i", "--input",   required=True, type=Path, help="Input TSV with accession/fetch_id/protein_id columns")
    ap.add_argument("-o", "--output",  required=True, type=Path, help="Output FASTA")
    ap.add_argument("--report",        required=True, type=Path, help="Output fetch report TSV")
    ap.add_argument("--sleep",         type=float, default=None, help="Minimum seconds between batch request starts")
    ap.add_argument("--workers",       type=int,   default=None, help="Concurrent batch workers; default: 4 with API key, 2 without")
    ap.add_argument("--batch-size",    type=int,   default=BATCH_SIZE, dest="batch_size",
                                       help=f"IDs per NCBI POST request (default: {BATCH_SIZE})")
    ap.add_argument("--email",         default="", help="NCBI email")
    ap.add_argument("--api_key",       default="", help="NCBI API key")
    ap.add_argument("--retries",       type=int,   default=3)
    args = ap.parse_args()

    rows = read_input(args.input)
    if not rows:
        die(f"no accessions to fetch in {args.input}")

    seen_protein_ids: set[str] = set()
    for row in rows:
        pid = row["protein_id"]
        if pid in seen_protein_ids:
            die(f"duplicate protein_id in exact fetch input: {pid}")
        seen_protein_ids.add(pid)

    sleep_sec = args.sleep if args.sleep is not None else (0.11 if args.api_key else 0.4)
    workers   = args.workers if args.workers is not None else (4 if args.api_key else 2)

    # Slice into batches; each batch is a list of (original_index, row) pairs.
    indexed = list(enumerate(rows, start=1))
    batch_size = max(1, args.batch_size)
    batches: list[list[tuple[int, Dict[str, str]]]] = [
        indexed[i : i + batch_size] for i in range(0, len(indexed), batch_size)
    ]
    workers = max(1, min(int(workers), len(batches)))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    info(f"Exact-fetch input rows: {len(rows)}")
    info(f"Output FASTA: {args.output}")
    info(f"Fetch report: {args.report}")
    info(f"Batch fetch: {len(batches)} batches x up to {batch_size} ids; workers={workers}; interval={sleep_sec}s")

    pacer   = RequestPacer(float(sleep_sec))
    results: Dict[int, Dict[str, object]] = {}
    completed = 0

    def run_batch(batch):
        return process_batch(batch, email=args.email or None, api_key=args.api_key or None,
                             retries=args.retries, pacer=pacer)

    try:
        if workers == 1:
            for batch in batches:
                for res in run_batch(batch):
                    results[res["idx"]] = res
                    completed += 1
                    if completed % 50 == 0 or completed == len(rows):
                        info(f"Fetched {completed}/{len(rows)} records")
        else:
            with futures.ThreadPoolExecutor(max_workers=workers) as ex:
                future_to_batch = {ex.submit(run_batch, b): b for b in batches}
                for fut in futures.as_completed(future_to_batch):
                    for res in fut.result():
                        results[res["idx"]] = res
                        completed += 1
                    if completed % 50 == 0 or completed == len(rows):
                        info(f"Fetched {completed}/{len(rows)} records")
    except Exception as exc:
        die(f"exact fetch failed: {exc}")

    fields = [
        "source_label", "protein_id", "fetch_id", "fetch_type", "source_accession",
        "returned_fasta_id", "returned_description", "sequence_length", "status", "note",
    ]
    with args.output.open("w", encoding="utf-8") as fout, \
         args.report.open("w", encoding="utf-8", newline="") as rout:
        writer = csv.DictWriter(rout, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for idx in range(1, len(rows) + 1):
            res = results[idx]
            fout.write(str(res["header"]) + "\n")
            for chunk in wrap(str(res["seq"])):
                fout.write(chunk + "\n")
            writer.writerow(res["report_row"])

    seq_count = sum(1 for line in args.output.open(encoding="utf-8") if line.startswith(">"))
    if seq_count != len(rows):
        die(f"output FASTA record count mismatch: expected {len(rows)}, observed {seq_count}")
    info(f"Exact fetch complete: {seq_count} sequences written")


if __name__ == "__main__":
    main()
