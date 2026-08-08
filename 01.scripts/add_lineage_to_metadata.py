#!/usr/bin/env python3

"""
add_lineage_to_metadata.py
──────────────────────────
maniFasta post-build step — taxonomic lineage enrichment.

Reads a per-protein metadata/manifest TSV (anything carrying a protein-id
column and an ncbi_taxid column), resolves each distinct taxid to a full
ranked NCBI lineage, and appends lineage columns to every row. The enriched
TSV is what plot_taxonomy / taxonomy_sunburst.html consumes.

Lineage source (two interchangeable paths):
  1. --lineage-cache CACHE.tsv
        Reuse a cache already produced by Fetch_TaxonomyLineage_INFO.sh
        (same column schema). No network needed. Preferred on the cluster.
  2. (no cache, or cache missing some taxids)
        Fetch the missing taxids directly from NCBI Taxonomy via
        epost + efetch, exactly like Fetch_TaxonomyLineage_INFO.sh. Requires
        outbound access to eutils.ncbi.nlm.nih.gov. Supply --email/--api-key
        to lift the rate limit. Fetched rows can be saved with --write-cache.

Appended columns (prefix configurable, default "tax_"):
  tax_own_rank tax_own_name tax_domain_or_realm tax_kingdom tax_phylum
  tax_class tax_order tax_family tax_genus tax_species tax_lineage_path

domain_or_realm merges cellular 'superkingdom' (Bacteria/Archaea/Eukaryota)
and viral 'realm' into one column, matching Fetch_TaxonomyLineage_INFO.sh.

Rows whose taxid is blank, non-numeric, or the NCBI "unidentified" default
(32644) are labelled as --unclassified-label so they still bucket cleanly in
the graphic instead of vanishing.

Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOLNAME = "maniFasta_addLineage"

# The NCBI "unidentified" taxid build_db.py assigns when nothing resolved.
UNIDENTIFIED_TAXID = "32644"

# Cache schema produced by Fetch_TaxonomyLineage_INFO.sh, in order.
CACHE_COLS = [
    "taxid", "own_rank", "own_name", "domain_or_realm", "kingdom", "phylum",
    "class", "order", "family", "genus", "species", "full_lineage_path",
]
# Ranks pulled from LineageEx (everything except domain_or_realm, handled separately).
RANK_COLS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def _ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def info(msg: str) -> None:
    print(f"[INFO] [{_ts()}] {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"[WARN] [{_ts()}] {msg}", file=sys.stderr)


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: [{_ts()}] {msg}", file=sys.stderr)
    raise SystemExit(code)


# ── TSV helpers ──────────────────────────────────────────────────────────────

def read_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    csv.field_size_limit(min(sys.maxsize, 2**31 - 1))
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    if not header:
        die(f"no header in {path}")
    return header, rows


def pick_col(header: List[str], *candidates: str) -> Optional[str]:
    low = {h.lower(): h for h in header}
    for c in candidates:
        if c.lower() in low:
            return low[c.lower()]
    return None


# ── Lineage cache ────────────────────────────────────────────────────────────

def load_cache(path: Path) -> Dict[str, Dict[str, str]]:
    """taxid -> {own_rank, own_name, domain_or_realm, <ranks>, full_lineage_path}."""
    cache: Dict[str, Dict[str, str]] = {}
    header, rows = read_tsv(path)
    missing = [c for c in ("taxid", "full_lineage_path") if c not in header]
    if missing:
        die(f"lineage cache {path} missing column(s): {missing}; "
            f"expected schema {CACHE_COLS}")
    for r in rows:
        tid = (r.get("taxid") or "").strip()
        if tid:
            cache[tid] = r
    info(f"Loaded {len(cache)} cached lineages from {path}")
    return cache


def write_cache(path: Path, cache: Dict[str, Dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CACHE_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for tid in sorted(cache, key=lambda x: int(x) if x.isdigit() else 0):
            row = cache[tid]
            w.writerow({c: row.get(c, "") for c in CACHE_COLS})
    info(f"Wrote {len(cache)} lineages to cache {path}")


# ── NCBI Taxonomy fetch (epost + efetch), parsed like the shell tool ─────────

def _http_post(url: str, params: Dict[str, str], *, retries: int, sleep_sec: float) -> str:
    data = urllib.parse.urlencode(params).encode("utf-8")
    last = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, method="POST",
                                         headers={"User-Agent": TOOLNAME})
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            warn(f"NCBI attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"NCBI request failed after {retries} attempts ({url}): {last}")


def parse_taxaset(xml_text: str) -> Dict[str, Dict[str, str]]:
    """Parse an efetch TaxaSet into {taxid: cache-row dict}."""
    out: Dict[str, Dict[str, str]] = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        warn(f"TaxaSet XML parse error: {e}")
        return out

    for taxon in root.findall("Taxon"):
        taxid = taxon.findtext("TaxId", "")
        name = taxon.findtext("ScientificName", "")
        own_rank = taxon.findtext("Rank", "no rank")

        rank_map: Dict[str, str] = {}
        path: List[str] = []
        lineage_ex = taxon.find("LineageEx")
        if lineage_ex is not None:
            for t in lineage_ex.findall("Taxon"):
                r = t.findtext("Rank", "no rank")
                n = t.findtext("ScientificName", "")
                if n:
                    path.append(n)
                if r in ("superkingdom", "domain", "realm"):
                    rank_map.setdefault("domain_or_realm", n)
                elif r in RANK_COLS:
                    rank_map[r] = n
        if name:
            path.append(name)
        if own_rank in ("superkingdom", "domain", "realm"):
            rank_map.setdefault("domain_or_realm", name)
        elif own_rank in RANK_COLS:
            rank_map[own_rank] = name

        row = {
            "taxid": taxid, "own_rank": own_rank, "own_name": name,
            "domain_or_realm": rank_map.get("domain_or_realm", ""),
            "full_lineage_path": "|".join(path),
        }
        for r in RANK_COLS:
            row[r] = rank_map.get(r, "")
        if taxid:
            out[taxid] = row
    return out


def fetch_lineages(taxids: List[str], *, email: str, api_key: str,
                   batch_size: int, sleep_sec: float, retries: int) -> Dict[str, Dict[str, str]]:
    resolved: Dict[str, Dict[str, str]] = {}
    common = {"tool": TOOLNAME}
    if email:
        common["email"] = email
    if api_key:
        common["api_key"] = api_key

    batches = [taxids[i:i + batch_size] for i in range(0, len(taxids), batch_size)]
    info(f"Fetching {len(taxids)} taxids from NCBI in {len(batches)} batch(es) "
         f"(size {batch_size}, interval {sleep_sec}s)")

    for n, batch in enumerate(batches, start=1):
        epost = _http_post(f"{EUTILS}/epost.fcgi",
                           {**common, "db": "taxonomy", "id": ",".join(batch)},
                           retries=retries, sleep_sec=sleep_sec)
        we = _between(epost, "<WebEnv>", "</WebEnv>")
        qk = _between(epost, "<QueryKey>", "</QueryKey>")
        if not (we and qk):
            warn(f"batch {n}: could not parse WebEnv/QueryKey; skipping {len(batch)} taxids")
            time.sleep(sleep_sec)
            continue
        time.sleep(sleep_sec)
        xml_text = _http_post(f"{EUTILS}/efetch.fcgi",
                              {**common, "db": "taxonomy", "WebEnv": we,
                               "query_key": qk, "retmode": "xml"},
                              retries=retries, sleep_sec=sleep_sec)
        got = parse_taxaset(xml_text)
        resolved.update(got)
        info(f"batch {n}/{len(batches)}: resolved {len(got)}/{len(batch)}")
        time.sleep(sleep_sec)
    return resolved


def _between(text: str, a: str, b: str) -> str:
    i = text.find(a)
    if i < 0:
        return ""
    j = text.find(b, i + len(a))
    if j < 0:
        return ""
    return text[i + len(a):j].strip()


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True, type=Path,
                    help="Input metadata/manifest TSV (per-protein rows).")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output enriched metadata TSV.")
    ap.add_argument("--lineage-cache", type=Path, default=None,
                    help="Existing lineage cache (Fetch_TaxonomyLineage_INFO.sh schema) to reuse.")
    ap.add_argument("--write-cache", type=Path, default=None,
                    help="Write the merged lineage cache (cached + freshly fetched) here.")
    ap.add_argument("--no-fetch", action="store_true",
                    help="Do not contact NCBI; only use --lineage-cache (unresolved taxids -> unclassified).")
    ap.add_argument("--id-col", default=None, help="Override the protein-id column name.")
    ap.add_argument("--taxid-col", default=None, help="Override the taxid column name.")
    ap.add_argument("--prefix", default="tax_", help="Prefix for appended lineage columns (default tax_).")
    ap.add_argument("--unclassified-label", default="Unclassified",
                    help="Label for rows with no resolvable lineage (default Unclassified).")
    ap.add_argument("--email", default="", help="NCBI email (raises rate limit).")
    ap.add_argument("--api-key", default="", dest="api_key", help="NCBI API key (raises rate limit).")
    ap.add_argument("--batch-size", type=int, default=300, dest="batch_size")
    ap.add_argument("--sleep", type=float, default=None,
                    help="Seconds between NCBI requests; default 0.11 with api-key else 0.4.")
    ap.add_argument("--retries", type=int, default=4)
    args = ap.parse_args()

    header, rows = read_tsv(args.metadata)

    id_col = args.id_col or pick_col(header, "protein_id", "Protein_ID", "id")
    taxid_col = args.taxid_col or pick_col(header, "ncbi_taxid", "taxid", "taxID", "TaxonID")
    if not taxid_col:
        die(f"could not find a taxid column in {args.metadata}; columns: {header}. "
            f"Use --taxid-col to set it.")
    if not id_col:
        warn(f"no protein-id column detected in {header}; continuing (id used only for logging).")
    info(f"Rows: {len(rows)} | id-col: {id_col or 'NA'} | taxid-col: {taxid_col}")

    # Distinct, resolvable taxids (numeric, not the unidentified default).
    def is_resolvable(tid: str) -> bool:
        return tid.isdigit() and tid != UNIDENTIFIED_TAXID

    distinct = sorted({(r.get(taxid_col) or "").strip() for r in rows
                       if is_resolvable((r.get(taxid_col) or "").strip())},
                      key=lambda x: int(x))
    n_unclassified_rows = sum(1 for r in rows
                              if not is_resolvable((r.get(taxid_col) or "").strip()))
    info(f"Distinct resolvable taxids: {len(distinct)} | "
         f"rows that will be {args.unclassified_label}: {n_unclassified_rows}")

    cache: Dict[str, Dict[str, str]] = {}
    if args.lineage_cache:
        cache = load_cache(args.lineage_cache)

    need = [t for t in distinct if t not in cache]
    if need and not args.no_fetch:
        sleep_sec = args.sleep if args.sleep is not None else (0.11 if args.api_key else 0.4)
        fetched = fetch_lineages(need, email=args.email, api_key=args.api_key,
                                 batch_size=max(1, args.batch_size),
                                 sleep_sec=sleep_sec, retries=args.retries)
        cache.update(fetched)
    elif need and args.no_fetch:
        warn(f"{len(need)} taxids absent from cache and --no-fetch set; "
             f"they will be {args.unclassified_label}.")

    if args.write_cache:
        write_cache(args.write_cache, cache)

    # Append columns.
    p = args.prefix
    new_cols = [p + c for c in
                ["own_rank", "own_name", "domain_or_realm", "kingdom", "phylum",
                 "class", "order", "family", "genus", "species", "lineage_path"]]
    out_header = header + [c for c in new_cols if c not in header]

    def unclassified_row() -> Dict[str, str]:
        lab = args.unclassified_label
        return {
            p + "own_rank": "", p + "own_name": lab, p + "domain_or_realm": lab,
            p + "kingdom": "", p + "phylum": "", p + "class": "", p + "order": "",
            p + "family": "", p + "genus": "", p + "species": "",
            p + "lineage_path": lab,
        }

    n_resolved = n_unresolved = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_header, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            tid = (r.get(taxid_col) or "").strip()
            lin = cache.get(tid) if is_resolvable(tid) else None
            if lin:
                r[p + "own_rank"] = lin.get("own_rank", "")
                r[p + "own_name"] = lin.get("own_name", "")
                r[p + "domain_or_realm"] = lin.get("domain_or_realm", "")
                for rk in RANK_COLS:
                    r[p + rk] = lin.get(rk, "")
                r[p + "lineage_path"] = lin.get("full_lineage_path", "")
                n_resolved += 1
            else:
                r.update(unclassified_row())
                n_unresolved += 1
            w.writerow(r)

    info(f"Enrichment complete: {n_resolved} rows with lineage, "
         f"{n_unresolved} {args.unclassified_label}.")
    info(f"Enriched metadata: {args.out}")


if __name__ == "__main__":
    main()
