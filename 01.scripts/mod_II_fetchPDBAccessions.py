#!/usr/bin/env python3

"""
mod_II_fetchPDBAccessions.py
─────────────────────
maniFasta Module 2 (accession fetch) — RCSB PDB backend.

Parallels mod_II_fetchUniProtAccessions.py but retrieves polymer sequences from the
RCSB PDB FASTA service:
    https://www.rcsb.org/fasta/entry/{ENTRY_ID}
Selection: reads the normalized source plan and processes rows where
    module == mod_II  AND  options contains fetch_source=pdb
Alternatively a flat accession list is read with -i/--input.

ENTITY EXPANSION (the one structural difference from the other fetchers):
  A PDB entry returns one FASTA record PER POLYMER ENTITY, so a single input
  accession can yield several output sequences. Two input forms are supported:

    entry-level   e.g. 4HHB        -> ALL polymer entities are emitted
                                       (4HHB_1, 4HHB_2, ...).
    entity-level  e.g. 4HHB_1      -> only that one entity is emitted.

  protein_id assignment keeps output ids unique:
    - entity-level input: protein_id = your protein_id if given, else the
      returned entity token (e.g. 4HHB_1).
    - entry-level input, N entities: if you supplied an explicit protein_id it
      is suffixed per entity (myid_1, myid_2, ...); otherwise the returned
      entity token is used. Single-entity entries keep the bare protein_id.

  Each request is fetched once per ENTRY (results cached) so multiple entity
  rows for the same entry do not re-hit the server.

Taxonomy:
  RCSB FASTA headers carry organism and taxid, e.g.
    >4HHB_1|Chains A, C|Hemoglobin subunit alpha|Homo sapiens (9606)
  so genus/species/ncbi_taxid are harvested directly. Chimeric entities list
  multiple organisms; the first is used and the record is noted as a chimera.

Outputs (orchestration appends these to user_fasta_file / user_metadata_file):
  --output    FASTA, first token rewritten to protein_id, provenance tokens added
  --metadata  user-metadata TSV (columns match build_db.py load_user_metadata())
  --report    fetch report TSV

With --allow-missing, accessions RCSB does not return are skipped (recorded in
the report with status=MISSING and optionally appended to --missing-out)
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
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

FASTA_BASE = "https://www.rcsb.org/fasta/entry"
TOOLNAME = "maniFasta_fetchPDBAccessions"

# Trailing _<digits> marks an entity-level id (4HHB_1, PDB_00001ABC_2); a bare
# entry id (4HHB, PDB_00001ABC) does not match and is treated as entry-level.
_ENTITY_RE = re.compile(r"^(?P<entry>.+?)_(?P<entity>\d+)$")
# Organism + NCBI taxid in a trailing "Name (12345)" fragment of the header.
_ORG_TAX_RE = re.compile(r"([^()|]+?)\s*\((\d+)\)")


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


def split_id(token: str) -> Tuple[str, Optional[str]]:
    """Return (entry_id, entity_number_or_None) for a PDB id token."""
    m = _ENTITY_RE.match(token)
    if m:
        return m.group("entry"), m.group("entity")
    return token, None


def parse_header_meta(header: str) -> Dict[str, str]:
    """Extract description, organism, taxid, chimera-flag from an RCSB header.

    Header shape: 'ENTITYID|Chains A, C|Molecule name|Organism (taxid)[, Org2 (tax2)]'
    """
    fields = header.split("|")
    description = fields[2].strip() if len(fields) >= 3 else clean_desc(header)
    orgs = _ORG_TAX_RE.findall(fields[-1] if fields else header)
    if orgs:
        organism, taxid = orgs[0][0].strip(), orgs[0][1].strip()
    else:
        organism, taxid = "", ""
    parts = organism.split(None, 1)
    return {
        "description": description or "NA",
        "organism": organism,
        "taxid": taxid,
        "genus": parts[0] if parts else "NA",
        "species": parts[1] if len(parts) > 1 else "NA",
        "chimera": "chimera_multi_organism" if len(orgs) > 1 else "",
    }


# ── HTTP ──────────────────────────────────────────────────────────────────────

def http_get(url: str, *, retries: int) -> str:
    headers = {"User-Agent": TOOLNAME}
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 404:
                return ""  # unknown entry; no record
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


def fetch_entry(entry: str, *, retries: int) -> str:
    url = f"{FASTA_BASE}/{urllib.parse.quote(entry)}"
    return http_get(url, retries=retries)


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header, chunks = None, []
    for raw in (text or "").splitlines():
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


# ── Plan / list reading (mirrors mod_II_fetchUniProtAccessions.py) ────────────────────

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


def _rows_to_requests(header: List[str], data: List[List[str]], default_label: str) -> List[Dict[str, str]]:
    i_acc = col_index(header, "accession")
    i_pid = col_index(header, "protein_id")
    i_fid = col_index(header, "fetch_id")
    i_sacc = col_index(header, "source_accession")
    i_label = col_index(header, "source_label")
    if i_acc < 0 and i_fid < 0:
        die(f"input needs an 'accession' or 'fetch_id' column. Found: {header}")
    out: List[Dict[str, str]] = []
    for row in data:
        def get(idx: int) -> str:
            return (row[idx].strip() if 0 <= idx < len(row) else "")
        acc = get(i_acc)
        fid = get(i_fid) or acc
        if not fid:
            continue
        pid = get(i_pid)
        out.append({
            "source_label": get(i_label) or default_label,
            # Empty protein_id signals "no explicit id" so expansion can use the
            # returned entity token. A non-empty value is treated as explicit.
            "protein_id": pid,
            "fetch_id": fid,
            "source_accession": get(i_sacc) or acc or fid,
        })
    return out


def collect_requests(plan_path: Path) -> List[Dict[str, str]]:
    """Gather fetch requests from all pdb mod_II rows."""
    out: List[Dict[str, str]] = []
    with plan_path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("enabled") or "").strip().lower() != "true":
                continue
            if (r.get("status") or "").strip().lower() == "disabled":
                continue
            if (r.get("module") or "").strip() != "mod_II":
                continue
            if parse_options(r.get("options") or "").get("fetch_source", "ncbi").lower() != "pdb":
                continue
            label = (r.get("source_label") or "PDB").strip()
            src = Path((r.get("input_list") or "").strip())
            if not src.is_file():
                die(f"input_list not found for '{label}': {src}")
            header, data = read_clean_table(src)
            out.extend(_rows_to_requests(header, data, label))
    return out


def collect_requests_from_list(list_path: Path) -> List[Dict[str, str]]:
    if not list_path.is_file():
        die(f"--input list not found: {list_path}")
    header, data = read_clean_table(list_path)
    return _rows_to_requests(header, data, "PDB")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", type=Path, help="Normalized source plan (selects pdb mod_II rows)")
    ap.add_argument("-i", "--input", type=Path, help="Flat accession list (alternative to --plan)")
    ap.add_argument("-o", "--output", required=True, type=Path, help="FASTA out")
    ap.add_argument("--metadata", required=True, type=Path, help="user-metadata TSV out")
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sleep", type=float, default=0.4)
    ap.add_argument("--retries", type=int, default=5)
    # accepted for CLI symmetry; PDB is fetched one entry per request.
    ap.add_argument("--batch-size", type=int, default=1, dest="batch_size")
    ap.add_argument("--api-key", default="")  # unused; RCSB FASTA needs no key
    ap.add_argument("--allow-missing", action="store_true",
                    help="Warn and skip accessions RCSB does not return instead of aborting.")
    ap.add_argument("--missing-out", type=Path, default=None,
                    help="Append skipped accessions to this TSV (source_label, accession, fetch_type, reason).")
    args = ap.parse_args()
    if bool(args.plan) == bool(args.input):
        die("provide exactly one of --plan or -i/--input")

    requests = collect_requests(args.plan) if args.plan else collect_requests_from_list(args.input)
    if not requests:
        info("No PDB accession rows in plan; nothing to do.")
        return

    info(f"PDB accession requests: {len(requests)}")

    # Fetch each distinct entry once; cache parsed entity records by uppercased
    # entity token (e.g. 4HHB_1) and group them under the entry.
    entry_cache: Dict[str, List[Tuple[str, str, str]]] = {}  # entry -> [(entity_token, hdr, seq)]
    entry_order: List[str] = []
    for r in requests:
        entry, _ent = split_id(r["fetch_id"])
        entry_u = entry.upper()
        if entry_u not in entry_cache:
            entry_order.append(entry_u)
            entry_cache[entry_u] = []  # placeholder so we only queue once

    for i, entry_u in enumerate(entry_order, start=1):
        text = fetch_entry(entry_u, retries=args.retries)
        recs: List[Tuple[str, str, str]] = []
        for hdr, seq in parse_fasta(text):
            tok = hdr.split("|")[0].strip() if "|" in hdr else (hdr.split()[0] if hdr.split() else "")
            recs.append((tok, hdr, seq))
        entry_cache[entry_u] = recs
        if i % 25 == 0 or i == len(entry_order):
            info(f"Fetched {i}/{len(entry_order)} entries")
        time.sleep(args.sleep)

    # Expand each request into one or more output records.
    out_records: List[Dict[str, object]] = []
    missing_reqs: List[Dict[str, str]] = []
    seen_pids: set[str] = set()

    def _unique(pid: str) -> str:
        base, n = pid, 2
        while pid in seen_pids:
            pid = f"{base}__{n}"
            n += 1
        seen_pids.add(pid)
        return pid

    for r in requests:
        entry, ent = split_id(r["fetch_id"])
        entry_u = entry.upper()
        recs = entry_cache.get(entry_u) or []
        explicit_pid = r["protein_id"].strip()

        if not recs:
            missing_reqs.append(r)
            continue

        if ent is not None:
            # entity-level: select the matching entity token (ENTRY_ENT).
            want = f"{entry_u}_{ent}"
            match = next((t for t in recs if t[0].upper() == want), None)
            if match is None:
                missing_reqs.append(r)
                continue
            chosen = [match]
            entry_level = False
        else:
            chosen = recs
            entry_level = True

        multi = len(chosen) > 1
        for (tok, hdr, seq) in chosen:
            hm = parse_header_meta(hdr)
            ent_num = (_ENTITY_RE.match(tok).group("entity")
                       if _ENTITY_RE.match(tok) else "1")
            if explicit_pid:
                pid = f"{explicit_pid}_{ent_num}" if (entry_level and multi) else explicit_pid
            else:
                pid = tok  # returned entity token, globally unique in PDB
            pid = _unique(clean_token(pid))

            out_records.append({
                "protein_id": pid,
                "source_label": r["source_label"],
                "fetch_id": r["fetch_id"],
                "source_accession": r["source_accession"],
                "returned_fasta_id": tok,
                "seq": seq,
                "meta": hm,
            })

    missing_ids = [r["fetch_id"] for r in missing_reqs]
    if missing_ids:
        if not args.allow_missing:
            die(f"{len(missing_ids)} accession(s) not returned by RCSB: {missing_ids[:25]}"
                + (" ..." if len(missing_ids) > 25 else "")
                + "\n  Re-run with --allow-missing to skip these and continue.")
        warn(f"{len(missing_ids)} accession(s) not returned by RCSB; skipping (--allow-missing): "
             f"{missing_ids[:25]}" + (" ..." if len(missing_ids) > 25 else ""))
        if args.missing_out:
            args.missing_out.parent.mkdir(parents=True, exist_ok=True)
            new = not args.missing_out.exists()
            with args.missing_out.open("a", encoding="utf-8", newline="") as mo:
                if new:
                    mo.write("source_label\taccession\tfetch_type\treason\n")
                for r in missing_reqs:
                    mo.write(f"{r['source_label']}\t{r['fetch_id']}\tpdb_entity\tnot_returned_by_rcsb\n")

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

        for rec in out_records:
            hm = rec["meta"]  # type: ignore[assignment]
            pid = rec["protein_id"]
            seq = rec["seq"]
            tok = rec["returned_fasta_id"]
            desc = hm["description"]

            out_header = (
                f">{clean_token(str(pid))} "
                f"fetch_id={clean_token(str(rec['fetch_id']))} "
                f"fetch_type=pdb_entity "
                f"source_accession={clean_token(str(rec['source_accession']))} "
                f"source_label={clean_token(str(rec['source_label']))} "
                f"returned_fasta_id={clean_token(str(tok))}"
            )
            if clean_desc(str(desc)):
                out_header += f" {clean_desc(str(desc))}"
            fout.write(out_header + "\n")
            for chunk in wrap(str(seq)):
                fout.write(chunk + "\n")

            mwr.writerow({
                "protein_id": clean_tsv(str(pid)),
                "source": clean_tsv(str(rec["source_label"])),
                "ncbi_taxid": hm["taxid"] or "",
                "genus": clean_tsv(hm["genus"]),
                "species": clean_tsv(hm["species"]),
                "strain_or_gene": "NA",
                "assembly": "NA",
                "description": clean_tsv(str(desc)),
                "fetch_id": clean_tsv(str(rec["fetch_id"])),
                "fetch_type": "pdb_entity",
                "source_accession": clean_tsv(str(rec["source_accession"])),
                "returned_fasta_id": clean_tsv(str(tok)),
            })
            rwr.writerow({
                "source_label": rec["source_label"], "protein_id": pid,
                "fetch_id": rec["fetch_id"], "fetch_type": "pdb_entity",
                "source_accession": rec["source_accession"], "returned_fasta_id": tok,
                "returned_description": desc, "sequence_length": len(str(seq)),
                "status": "OK", "note": hm["chimera"],
            })

        for r in missing_reqs:
            rwr.writerow({
                "source_label": r["source_label"], "protein_id": r["protein_id"] or r["fetch_id"],
                "fetch_id": r["fetch_id"], "fetch_type": "pdb_entity",
                "source_accession": r["source_accession"], "returned_fasta_id": "NA",
                "returned_description": "", "sequence_length": 0,
                "status": "MISSING", "note": "not_returned_by_rcsb",
            })

    n = sum(1 for ln in args.output.open(encoding="utf-8") if ln.startswith(">"))
    if n != len(out_records):
        die(f"output record count mismatch: expected {len(out_records)}, wrote {n}")
    info(f"PDB accession fetch complete: {n} sequences from "
         f"{len(requests) - len(missing_reqs)} resolved request(s) "
         f"({len(missing_reqs)} missing)")
    info(f"FASTA: {args.output}")
    info(f"Metadata: {args.metadata}")
    info(f"Report: {args.report}")


if __name__ == "__main__":
    main()
