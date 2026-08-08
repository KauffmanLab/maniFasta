#!/usr/bin/env python3

"""
mod_II_fetchUniProtAccessions.py
─────────────────────────
maniFasta Module 2 (accession fetch) — UniProt backend.

Parallels mod_II_fetchProteinAccessions_batch.py but retrieves UniProtKB entries by
accession instead of NCBI protein records. Selection: reads the normalized
source plan and processes rows where
    module == mod_II  AND  options contains fetch_source=uniprot

Per input_list row the fetch_id is a UniProt accession (e.g. P02666, P02666-2).
Returned records are matched back to requests by the accession token in the
returned header (db|ACC|ENTRY), with version/isoform-stripped fallback — so
response ordering is irrelevant.

Outputs (the orchestration appends these to user_fasta_file / user_metadata_file,
exactly like the NCBI protein path):

  --output    FASTA, first token rewritten to protein_id, provenance tokens added
  --metadata  user-metadata TSV derived from UniProt headers (OS/OX/GN)
              columns match build_db.py load_user_metadata()
  --report    fetch report TSV

With --allow-missing, accessions UniProt does not return are skipped (recorded
in the report with status=MISSING and optionally appended to --missing-out)
instead of aborting the build.

Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REST = "https://rest.uniprot.org"
TOOLNAME = "maniFasta_fetchUniProtAccessions"
BATCH_SIZE = 100  # accessions per UniProt request


def _ts() -> str:
    return datetime.datetime.now().strftime("%H:%M:%S")


def info(msg: str) -> None:
    print(f"[INFO] [{_ts()}] {msg}", file=sys.stderr)


def warn(msg: str) -> None:
    print(f"[WARN] [{_ts()}] {msg}", file=sys.stderr)


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: [{_ts()}] {msg}", file=sys.stderr)
    raise SystemExit(code)


def parse_options(raw: str) -> Dict[str, str]:
    opts: Dict[str, str] = {}
    for item in (raw or "").split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        k, v = item.split("=", 1)
        opts[k.strip().lower()] = v.strip()
    return opts


def clean_token(value: str) -> str:
    value = (value or "NA").strip() or "NA"
    return "_".join(value.replace("\r", " ").replace("\n", " ").split())


def clean_desc(value: str) -> str:
    return " ".join((value or "").replace("\r", " ").replace("\n", " ").split())


def clean_tsv(value: str) -> str:
    value = str(value if value is not None else "").replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()
    return value or "NA"


def wrap(seq: str, width: int = 60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


# ── HTTP ──────────────────────────────────────────────────────────────────────

def http_get(url: str, *, retries: int, api_key: Optional[str], timeout: int = 180) -> str:
    headers = {"User-Agent": TOOLNAME}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 429:
                wait = float(e.headers.get("Retry-After", 2 ** attempt))
                warn(f"429 rate limited; sleeping {wait}s")
                time.sleep(wait)
                continue
            warn(f"HTTP {e.code} (attempt {attempt}/{retries})")
        except Exception as e:  # noqa: BLE001
            last_err = e
            warn(f"GET attempt {attempt}/{retries} failed: {e}")
        if attempt < retries:
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"GET failed after {retries} attempts: {last_err}")


def fetch_accessions(accs: List[str], *, retries: int, api_key: Optional[str]) -> str:
    # UniProt's REST API takes bare accessions (no NCBI-style ".N" version
    # suffix); a versioned token like "P00785.4" gets the whole batch
    # rejected with HTTP 400. Strip it for the query only -- the unstripped
    # fetch_id is still what we match results back against below.
    query_accs = [a.split(".")[0] for a in accs]
    params = {"accessions": ",".join(query_accs), "format": "fasta"}
    url = f"{REST}/uniprotkb/accessions?{urllib.parse.urlencode(params)}"
    return http_get(url, retries=retries, api_key=api_key)


# ── Taxon-scope bulk fetch (uniprot_scope=taxon) ──────────────────────────────
# A taxon row pulls every UniProtKB protein under a taxid in one stream, instead
# of one-accession-at-a-time. Reviewed (Swiss-Prot) only by default; opt into
# the full set (incl. TrEMBL) with uniprot_include_trembl=true. Each returned
# protein keeps its OWN organism/taxid/gene from its header, so a multi-species
# taxon (genus/family) retains correct per-protein taxonomy downstream. There is
# deliberately no result-size cap.

def build_taxon_query(taxid: str, *, include_trembl: bool, include_isoform: bool) -> Tuple[str, str]:
    """Return (stream_url, human_readable_query) for a taxid bulk fetch."""
    query = f"(taxonomy_id:{taxid})"
    if not include_trembl:
        query += " AND (reviewed:true)"
    params = {
        "query": query,
        "format": "fasta",
        "includeIsoform": "true" if include_isoform else "false",
    }
    return f"{REST}/uniprotkb/stream?{urllib.parse.urlencode(params)}", query


def fetch_taxon_stream(taxid: str, *, include_trembl: bool, include_isoform: bool,
                       retries: int, api_key: Optional[str]) -> Tuple[List[Tuple[str, str]], str]:
    """Stream all UniProtKB proteins for a taxid; returns (records, query).

    Uses a long timeout because an unbounded taxonomy_id stream (especially with
    TrEMBL) can be large.
    """
    url, query = build_taxon_query(taxid, include_trembl=include_trembl,
                                   include_isoform=include_isoform)
    text = http_get(url, retries=retries, api_key=api_key, timeout=600)
    return parse_fasta(text), query


# ── FASTA / header parsing ────────────────────────────────────────────────────

def parse_fasta(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header, chunks = None, []
    for raw in text.splitlines():
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def header_accession(token: str) -> str:
    """db|ACC|ENTRY -> ACC ; otherwise the token itself."""
    fields = token.split("|")
    return fields[1] if len(fields) >= 3 else token


_RE_OS = re.compile(r"\bOS=(.+?)(?:\s+[A-Z]{2}=|$)")
_RE_OX = re.compile(r"\bOX=(\d+)")
_RE_GN = re.compile(r"\bGN=([^\s]+)")


def header_metadata(desc: str) -> dict:
    os_m = _RE_OS.search(desc)
    organism = os_m.group(1).strip() if os_m else ""
    parts = organism.split(None, 1)
    return {
        "taxid": (_RE_OX.search(desc).group(1) if _RE_OX.search(desc) else ""),
        "genus": parts[0] if parts else "NA",
        "species": parts[1] if len(parts) > 1 else "NA",
        "gene": (_RE_GN.search(desc).group(1) if _RE_GN.search(desc) else "NA"),
        "description": clean_desc(desc) or "NA",
    }


# ── Plan reading ──────────────────────────────────────────────────────────────

def read_clean_table(path: Path) -> tuple[List[str], List[List[str]]]:
    lines = [ln.rstrip("\n").replace("\r", "")
             for ln in path.open(encoding="utf-8", errors="replace")
             if ln.strip() and not ln.startswith("#")]
    if not lines:
        die(f"no data rows in {path}")
    header = [h.strip() for h in lines[0].split("\t")]
    return header, [ln.split("\t") for ln in lines[1:]]


def col_index(header: List[str], *names: str) -> int:
    low = [h.lower() for h in header]
    for n in names:
        if n.lower() in low:
            return low.index(n.lower())
    return -1


def collect_requests(plan_path: Path) -> List[Dict[str, str]]:
    """Gather fetch requests from all uniprot mod_II rows."""
    out: List[Dict[str, str]] = []
    with plan_path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("enabled") or "").strip().lower() != "true":
                continue
            if (r.get("status") or "").strip().lower() == "disabled":
                continue
            if (r.get("module") or "").strip() != "mod_II":
                continue
            _o = parse_options(r.get("options") or "")
            if _o.get("fetch_source", "ncbi").lower() != "uniprot":
                continue
            if _o.get("uniprot_scope", "").lower() == "taxon":
                continue  # taxon-scope rows are handled by collect_taxon_sources()

            label = (r.get("source_label") or "UNIPROT").strip()
            src = Path((r.get("input_list") or "").strip())
            if not src.is_file():
                die(f"input_list not found for '{label}': {src}")
            header, data = read_clean_table(src)
            i_acc = col_index(header, "accession")
            i_pid = col_index(header, "protein_id")
            i_fid = col_index(header, "fetch_id")
            i_sacc = col_index(header, "source_accession")
            if i_acc < 0 and i_fid < 0:
                die(f"'{label}': input_list needs an 'accession' or 'fetch_id' column. Found: {header}")
            for row in data:
                def get(idx: int) -> str:
                    return (row[idx].strip() if 0 <= idx < len(row) else "")
                acc = get(i_acc)
                fid = get(i_fid) or acc
                if not fid:
                    continue
                out.append({
                    "source_label": label,
                    "protein_id": get(i_pid) or fid,
                    "fetch_id": fid,
                    "source_accession": get(i_sacc) or acc or fid,
                })
    return out


def collect_requests_from_list(list_path: Path) -> List[Dict[str, str]]:
    """Gather requests from a flat accession list (as written by the staging
    partitioner). Same request shape as collect_requests()."""
    if not list_path.is_file():
        die(f"--input list not found: {list_path}")
    header, data = read_clean_table(list_path)
    i_acc   = col_index(header, "accession")
    i_pid   = col_index(header, "protein_id")
    i_fid   = col_index(header, "fetch_id")
    i_sacc  = col_index(header, "source_accession")
    i_label = col_index(header, "source_label")
    if i_acc < 0 and i_fid < 0:
        die(f"--input list needs an 'accession' or 'fetch_id' column. Found: {header}")
    out: List[Dict[str, str]] = []
    for row in data:
        def get(idx: int) -> str:
            return (row[idx].strip() if 0 <= idx < len(row) else "")
        acc = get(i_acc)
        fid = get(i_fid) or acc
        if not fid:
            continue
        out.append({
            "source_label": get(i_label) or "UNIPROT",
            "protein_id": get(i_pid) or fid,
            "fetch_id": fid,
            "source_accession": get(i_sacc) or acc or fid,
        })
    return out


def _truthy(v: str) -> bool:
    return (v or "").strip().lower() in ("true", "1", "yes", "on")


def collect_taxon_sources(plan_path: Path, *, default_trembl: bool = False,
                          default_isoform: bool = False) -> List[Dict[str, object]]:
    """Gather taxon-scope bulk sources from the plan: mod_II rows with
    fetch_source=uniprot AND uniprot_scope=taxon. Each such row's input_list
    supplies one or more taxids (a 'taxid' column); each taxid becomes one
    bulk stream.

    Reviewed/TrEMBL and isoform inclusion follow the global defaults
    (default_trembl / default_isoform, set from config via the CLI) UNLESS the
    row's options override them with uniprot_include_trembl= / uniprot_include_isoform=.
    """
    out: List[Dict[str, object]] = []
    with plan_path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("enabled") or "").strip().lower() != "true":
                continue
            if (r.get("status") or "").strip().lower() == "disabled":
                continue
            if (r.get("module") or "").strip() != "mod_II":
                continue
            opts = parse_options(r.get("options") or "")
            if opts.get("fetch_source", "ncbi").lower() != "uniprot":
                continue
            if opts.get("uniprot_scope", "").lower() != "taxon":
                continue

            label = (r.get("source_label") or "UNIPROT").strip()
            src = Path((r.get("input_list") or "").strip())
            if not src.is_file():
                die(f"input_list not found for '{label}': {src}")
            header, data = read_clean_table(src)
            i_tax = col_index(header, "taxid", "taxonID", "taxon_id", "ncbi_taxid")
            if i_tax < 0:
                die(f"'{label}': taxon-scope input_list needs a 'taxid' column. "
                    f"Found: {header}")
            # Per-row option overrides the global (config) default when present.
            include_trembl = (_truthy(opts["uniprot_include_trembl"])
                              if "uniprot_include_trembl" in opts else default_trembl)
            include_isoform = (_truthy(opts["uniprot_include_isoform"])
                               if "uniprot_include_isoform" in opts else default_isoform)
            for row in data:
                taxid = (row[i_tax].strip() if 0 <= i_tax < len(row) else "").replace(" ", "")
                if not taxid:
                    continue
                if not taxid.isdigit():
                    die(f"'{label}': non-numeric taxid {taxid!r} in {src}")
                out.append({
                    "source_label": label,
                    "taxid": taxid,
                    "include_trembl": include_trembl,
                    "include_isoform": include_isoform,
                })
    return out


# ── Main ────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", type=Path, help="Normalized source plan (selects uniprot mod_II rows)")
    ap.add_argument("-i", "--input", type=Path, help="Flat accession list (alternative to --plan)")
    ap.add_argument("-o", "--output", required=True, type=Path, help="FASTA out")
    ap.add_argument("--metadata", required=True, type=Path, help="user-metadata TSV out")
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sleep", type=float, default=0.4)
    ap.add_argument("--retries", type=int, default=5)
    ap.add_argument("--batch-size", type=int, default=BATCH_SIZE, dest="batch_size")
    ap.add_argument("--api-key", default="")
    ap.add_argument("--include-trembl", default="false", dest="include_trembl",
                    help="(taxon scope) global default: include TrEMBL, not just "
                         "reviewed Swiss-Prot. Set from config; a per-row "
                         "uniprot_include_trembl= option overrides it.")
    ap.add_argument("--include-isoform", default="false", dest="include_isoform",
                    help="(taxon scope) global default: include isoforms. Set from "
                         "config; a per-row uniprot_include_isoform= option overrides it.")
    ap.add_argument("--allow-missing", action="store_true",
                    help="Warn and skip accessions UniProt does not return instead of aborting.")
    ap.add_argument("--missing-out", type=Path, default=None,
                    help="Append skipped accessions to this TSV (source_label, accession, fetch_type, reason).")
    args = ap.parse_args()
    if bool(args.plan) == bool(args.input):
        die("provide exactly one of --plan or -i/--input")

    api_key = args.api_key or None
    requests = collect_requests(args.plan) if args.plan else collect_requests_from_list(args.input)
    taxon_sources = (collect_taxon_sources(
                         args.plan,
                         default_trembl=_truthy(args.include_trembl),
                         default_isoform=_truthy(args.include_isoform))
                     if args.plan else [])
    if not requests and not taxon_sources:
        info("No UniProt accession or taxon rows in plan; nothing to do.")
        return

    seen = set()
    for r in requests:
        if r["protein_id"] in seen:
            die(f"duplicate protein_id across uniprot accession inputs: {r['protein_id']}")
        seen.add(r["protein_id"])

    if requests:
        info(f"UniProt accession requests: {len(requests)}")
    batch_size = max(1, args.batch_size)

    # fetch_id -> (returned_header, seq), indexed by accession token + bare form
    returned: Dict[str, Tuple[str, str]] = {}
    for i in range(0, len(requests), batch_size):
        batch = requests[i:i + batch_size]
        text = fetch_accessions([r["fetch_id"] for r in batch], retries=args.retries, api_key=api_key)
        for hdr, seq in parse_fasta(text):
            tok = hdr.split()[0] if hdr.split() else ""
            acc = header_accession(tok)
            if acc:
                returned[acc] = (hdr, seq)
                returned.setdefault(acc.split("-")[0], (hdr, seq))  # isoform-bare
        info(f"Fetched {min(i + batch_size, len(requests))}/{len(requests)}")
        time.sleep(args.sleep)

    # A request is retrievable if its fetch_id (exact, isoform-bare, or
    # version-bare) is present in the returned index. Mirrors the lookup used
    # in the writing loop below so any kept request is guaranteed to resolve.
    def _is_missing(r: Dict[str, str]) -> bool:
        fid = r["fetch_id"]
        return (fid not in returned
                and fid.split("-")[0] not in returned
                and fid.split(".")[0] not in returned)

    missing_reqs = [r for r in requests if _is_missing(r)]
    kept_requests = [r for r in requests if not _is_missing(r)]
    missing_ids = [r["fetch_id"] for r in missing_reqs]

    if missing_ids:
        if not args.allow_missing:
            die(f"{len(missing_ids)} accession(s) not returned by UniProt: {missing_ids[:25]}"
                + (" ..." if len(missing_ids) > 25 else "")
                + "\n  Re-run with --allow-missing to skip these and continue.")
        warn(f"{len(missing_ids)} accession(s) not returned by UniProt; skipping (--allow-missing): "
             f"{missing_ids[:25]}" + (" ..." if len(missing_ids) > 25 else ""))
        if args.missing_out:
            args.missing_out.parent.mkdir(parents=True, exist_ok=True)
            new = not args.missing_out.exists()
            with args.missing_out.open("a", encoding="utf-8", newline="") as mo:
                if new:
                    mo.write("source_label\taccession\tfetch_type\treason\n")
                for r in missing_reqs:
                    mo.write(f"{r['source_label']}\t{r['fetch_id']}\tuniprot_accession\tnot_returned_by_uniprot\n")

    # ── Taxon-scope bulk streams ─────────────────────────────────────────────
    # Each taxon source pulls all UniProtKB proteins under its taxid. Every
    # returned protein keeps its own header-derived taxonomy. Dedup by accession
    # against accession-fetched protein_ids and across taxon streams.
    taxon_records: List[Dict[str, object]] = []
    for tsrc in taxon_sources:
        scope_desc = ("reviewed+TrEMBL" if tsrc["include_trembl"] else "reviewed only")
        if tsrc["include_isoform"]:
            scope_desc += ", +isoforms"
        info(f"[{tsrc['source_label']}] taxon fetch: taxonomy_id:{tsrc['taxid']} ({scope_desc})")
        recs, query = fetch_taxon_stream(
            str(tsrc["taxid"]),
            include_trembl=bool(tsrc["include_trembl"]),
            include_isoform=bool(tsrc["include_isoform"]),
            retries=args.retries, api_key=api_key)
        kept = 0
        dup = 0
        for hdr, seq in recs:
            tok = hdr.split()[0] if hdr.split() else ""
            acc = header_accession(tok)
            if not acc:
                continue
            if acc in seen:
                dup += 1
                continue
            seen.add(acc)
            parts = hdr.split()
            returned_id = parts[0] if parts else "NA"
            returned_desc = " ".join(parts[1:]) if len(parts) > 1 else ""
            taxon_records.append({
                "source_label": tsrc["source_label"],
                "protein_id": acc,
                "returned_id": returned_id,
                "returned_desc": returned_desc,
                "seq": seq,
                "query": query,
                "hm": header_metadata(returned_desc),
            })
            kept += 1
        info(f"  taxonomy_id:{tsrc['taxid']} -> {kept} sequences"
             + (f" ({dup} duplicate protein_id skipped)" if dup else ""))
        time.sleep(args.sleep)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    meta_cols = ["protein_id", "source", "ncbi_taxid", "genus", "species",
                 "strain_or_gene", "assembly", "description",
                 "fetch_id", "fetch_type", "source_accession", "returned_fasta_id"]
    rep_cols = ["source_label", "protein_id", "fetch_id", "fetch_type", "source_accession",
                "returned_fasta_id", "returned_description", "sequence_length", "status", "note"]

    with args.output.open("w", encoding="utf-8") as fout, \
         args.metadata.open("w", encoding="utf-8", newline="") as mout, \
         args.report.open("w", encoding="utf-8", newline="") as rout:

        mwr = csv.DictWriter(mout, fieldnames=meta_cols, delimiter="\t", lineterminator="\n")
        rwr = csv.DictWriter(rout, fieldnames=rep_cols, delimiter="\t", lineterminator="\n")
        mwr.writeheader()
        rwr.writeheader()

        for r in kept_requests:
            fid = r["fetch_id"]
            hdr, seq = returned.get(fid) or returned.get(fid.split("-")[0]) or returned[fid.split(".")[0]]
            parts = hdr.split()
            returned_id = parts[0] if parts else "NA"
            returned_desc = " ".join(parts[1:]) if len(parts) > 1 else ""
            hm = header_metadata(returned_desc)

            out_header = (
                f">{clean_token(r['protein_id'])} "
                f"fetch_id={clean_token(fid)} "
                f"fetch_type=uniprot_accession "
                f"source_accession={clean_token(r['source_accession'])} "
                f"returned_fasta_id={clean_token(returned_id)}"
            )
            if clean_desc(returned_desc):
                out_header += f" {clean_desc(returned_desc)}"
            fout.write(out_header + "\n")
            for chunk in wrap(seq):
                fout.write(chunk + "\n")

            mwr.writerow({
                "protein_id": clean_tsv(r["protein_id"]),
                "source": clean_tsv(r["source_label"]),
                "ncbi_taxid": hm["taxid"],
                "genus": clean_tsv(hm["genus"]),
                "species": clean_tsv(hm["species"]),
                "strain_or_gene": clean_tsv(hm["gene"]),
                "assembly": "NA",
                "description": clean_tsv(hm["description"]),
                "fetch_id": clean_tsv(fid),
                "fetch_type": "uniprot_accession",
                "source_accession": clean_tsv(r["source_accession"]),
                "returned_fasta_id": clean_tsv(returned_id),
            })
            rwr.writerow({
                "source_label": r["source_label"], "protein_id": r["protein_id"],
                "fetch_id": fid, "fetch_type": "uniprot_accession",
                "source_accession": r["source_accession"], "returned_fasta_id": returned_id,
                "returned_description": returned_desc, "sequence_length": len(seq),
                "status": "OK", "note": "",
            })

        # Record skipped accessions in the fetch report (FASTA/metadata omit them).
        for r in missing_reqs:
            rwr.writerow({
                "source_label": r["source_label"], "protein_id": r["protein_id"],
                "fetch_id": r["fetch_id"], "fetch_type": "uniprot_accession",
                "source_accession": r["source_accession"], "returned_fasta_id": "NA",
                "returned_description": "", "sequence_length": 0,
                "status": "MISSING", "note": "not_returned_by_uniprot",
            })

        # Taxon-scope records: same output shape, fetch_type=uniprot_taxon.
        # protein_id / source_accession are the discovered accession; the taxon
        # query is recorded in the report note for provenance.
        for tr in taxon_records:
            acc = str(tr["protein_id"])
            returned_id = str(tr["returned_id"])
            returned_desc = str(tr["returned_desc"])
            hm = tr["hm"]
            out_header = (
                f">{clean_token(acc)} "
                f"fetch_id={clean_token(acc)} "
                f"fetch_type=uniprot_taxon "
                f"source_accession={clean_token(acc)} "
                f"returned_fasta_id={clean_token(returned_id)}"
            )
            if clean_desc(returned_desc):
                out_header += f" {clean_desc(returned_desc)}"
            fout.write(out_header + "\n")
            for chunk in wrap(str(tr["seq"])):
                fout.write(chunk + "\n")

            mwr.writerow({
                "protein_id": clean_tsv(acc),
                "source": clean_tsv(str(tr["source_label"])),
                "ncbi_taxid": hm["taxid"],
                "genus": clean_tsv(hm["genus"]),
                "species": clean_tsv(hm["species"]),
                "strain_or_gene": clean_tsv(hm["gene"]),
                "assembly": "NA",
                "description": clean_tsv(hm["description"]),
                "fetch_id": clean_tsv(acc),
                "fetch_type": "uniprot_taxon",
                "source_accession": clean_tsv(acc),
                "returned_fasta_id": clean_tsv(returned_id),
            })
            rwr.writerow({
                "source_label": str(tr["source_label"]), "protein_id": acc,
                "fetch_id": acc, "fetch_type": "uniprot_taxon",
                "source_accession": acc, "returned_fasta_id": returned_id,
                "returned_description": returned_desc, "sequence_length": len(str(tr["seq"])),
                "status": "OK", "note": f"taxon_query={tr['query']}",
            })

    n = sum(1 for ln in args.output.open(encoding="utf-8") if ln.startswith(">"))
    expected = len(kept_requests) + len(taxon_records)
    if n != expected:
        die(f"output record count mismatch: expected {expected}, wrote {n}")
    info(f"UniProt fetch complete: {n} sequences "
         f"({len(kept_requests)} by accession, {len(taxon_records)} by taxon; "
         f"{len(missing_ids)} accession(s) missing)")
    info(f"FASTA: {args.output}")
    info(f"Metadata: {args.metadata}")
    info(f"Report: {args.report}")


if __name__ == "__main__":
    main()
