#!/usr/bin/env python3

"""
mod_II_fetchUniParcAccessions.py
─────────────────────────
maniFasta Module 2 (accession fetch) — UniParc backend.

Parallels mod_II_fetchUniProtAccessions.py but retrieves UniParc (UniProt Archive)
records by UPI accession (e.g. UPI0000000001). Selection: reads the normalized
source plan and processes rows where
    module == mod_II  AND  options contains fetch_source=uniparc
Alternatively a flat accession list is read with -i/--input.

NO METADATA FILE REQUIRED
─────────────────────────
Like the NCBI protein path (which auto-stubs metadata from the returned NCBI
title), this fetcher builds full taxonomy itself, so a bare UPI list needs no
metadata file. It does this by fetching the UniParc *entry JSON* rather than
the FASTA: the JSON carries the sequence AND the database cross-references,
each of which has an organism (taxonId + scientificName), proteinName, and
geneName. The bare UniParc FASTA header ('>UPI... status=active') has none of
this, which is why JSON is used.

  Endpoint: https://rest.uniprot.org/uniparc/{upi}.json

REPRESENTATIVE ORGANISM (inherent UniParc ambiguity)
────────────────────────────────────────────────────
A UPI is one unique sequence that may be seen in many organisms, so there is
no single "true" taxon. A representative is chosen deterministically
(--organism-pick, default 'best'):

  best   : restrict to the most authoritative active cross-references
           (reviewed Swiss-Prot > TrEMBL/UniProtKB > RefSeq > any), then the
           most frequent organism within that tier.
  common : the most frequent organism across all cross-references.
  first  : the organism on the first cross-reference (JSON order).

The number of distinct organisms is recorded in the report 'note' column
(n_organisms=K), flagged 'ambiguous_organism' when K>1, so the choice is never
silent. ncbi_taxid/genus/species reflect the chosen representative; supply a
metadata_file in the registry to override per-protein if you need exact taxa.
If no cross-reference carries an organism, commonTaxons is used as a last
resort; failing that, taxonomy is left blank (build_db defaults taxid to
32644 = unidentified).

Outputs (orchestration appends these to user_fasta_file / user_metadata_file):
  --output    FASTA, first token rewritten to protein_id, provenance tokens added
  --metadata  user-metadata TSV (columns match build_db.py load_user_metadata())
  --report    fetch report TSV

With --allow-missing, accessions UniParc does not return are skipped (recorded
in the report with status=MISSING and optionally appended to --missing-out)
instead of aborting the build.

Stdlib only.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

REST = "https://rest.uniprot.org"
TOOLNAME = "maniFasta_fetchUniParcAccessions"
UPI_RE = re.compile(r"^UPI[0-9A-Fa-f]+$")


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

def http_get(url: str, *, retries: int, api_key: Optional[str]) -> str:
    headers = {"User-Agent": TOOLNAME}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 404:
                return ""  # unknown UPI; no record
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


def fetch_entry_json(upi: str, *, retries: int, api_key: Optional[str]) -> Optional[dict]:
    """Fetch a UniParc entry as JSON. None if not found / unparseable."""
    url = f"{REST}/uniparc/{urllib.parse.quote(upi)}.json"
    text = http_get(url, retries=retries, api_key=api_key)
    if not text.strip():
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError as e:
        warn(f"JSON parse error for {upi}: {e}")
        return None


# ── Metadata derivation from UniParc entry JSON ───────────────────────────────

def _db_rank(database: str) -> int:
    d = (database or "").lower()
    if "swiss-prot" in d or "swissprot" in d:
        return 3
    if "trembl" in d:
        return 2
    if d.startswith("uniprotkb") or d == "uniprot":
        return 2
    if "refseq" in d:
        return 1
    return 0


def _org_fields(org: dict) -> tuple[str, str]:
    """Return (taxid, scientificName) from an organism object, tolerant of keys."""
    if not isinstance(org, dict):
        return "", ""
    taxid = org.get("taxonId", org.get("taxonomyId", org.get("taxId", "")))
    taxid = str(taxid).strip() if taxid not in (None, "") else ""
    name = (org.get("scientificName") or org.get("name") or "").strip()
    return taxid, name


def derive_metadata(entry: dict, *, pick: str) -> dict:
    """Derive taxonomy/description from a UniParc entry JSON.

    Returns a dict with: seq, taxid, genus, species, gene, description,
    organism, picked_db, n_organisms, note.
    """
    seq = ((entry.get("sequence") or {}).get("value") or "").strip()
    upi = (entry.get("uniParcId") or "").strip()

    xrefs = entry.get("uniParcCrossReferences") or entry.get("dbReference") or []
    # Candidate cross-references that carry an organism.
    cands = []  # (active, db_rank, taxid, sciname, proteinName, geneName)
    for x in xrefs:
        if not isinstance(x, dict):
            continue
        taxid, sciname = _org_fields(x.get("organism") or {})
        if not taxid and not sciname:
            continue
        active = bool(x.get("active", True))
        rank = _db_rank(x.get("database") or x.get("databaseType") or "")
        pname = (x.get("proteinName") or "").strip()
        gname = (x.get("geneName") or "").strip()
        cands.append((active, rank, taxid, sciname, pname, gname))

    note_bits: List[str] = []
    n_orgs = len({c[2] for c in cands if c[2]})

    chosen = None
    if cands:
        if pick == "first":
            chosen = cands[0]
        else:
            if pick == "best":
                top_key = max((c[0], c[1]) for c in cands)
                pool = [c for c in cands if (c[0], c[1]) == top_key]
            else:  # common
                pool = cands
            # Most frequent organism (by taxid) within the pool.
            tax_counts = Counter(c[2] for c in pool if c[2])
            if tax_counts:
                best_tax = tax_counts.most_common(1)[0][0]
                # representative xref for naming = highest db_rank for that taxon
                org_pool = [c for c in pool if c[2] == best_tax]
                chosen = max(org_pool, key=lambda c: (c[0], c[1]))
            else:
                chosen = max(pool, key=lambda c: (c[0], c[1]))
        note_bits.append(f"picked_db_rank={chosen[1]}")
        note_bits.append(f"n_organisms={n_orgs}")
        if n_orgs > 1:
            note_bits.append("ambiguous_organism")
    else:
        # Fallback: commonTaxons LCA, if present (can be broad).
        for ct in (entry.get("commonTaxons") or []):
            tid = str(ct.get("commonTaxonId") or "").strip()
            nm = (ct.get("commonTaxon") or "").strip()
            if tid or nm:
                chosen = (True, 0, tid, nm, "", "")
                note_bits.append("organism_from_commonTaxon")
                break
        if chosen is None:
            note_bits.append("no_organism_in_xrefs")

    if chosen is not None:
        _, _, taxid, sciname, pname, gname = chosen
    else:
        taxid, sciname, pname, gname = "", "", "", ""

    parts = sciname.split(None, 1)
    # build_db hard-fails any user row that has a non-default taxid but genus
    # or species == NA. A UniParc representative can legitimately be genus-level
    # or a commonTaxon LCA (e.g. "Bacteria"), which has no species. In that case
    # we must NOT assert a partial taxid in the authoritative columns, or the
    # whole build aborts. So we surface taxid + genus + species ONLY for a full
    # binomial; anything coarser is left unidentified (taxid blank -> build_db's
    # 32644 default) with the resolved value preserved in the report note for
    # auditing or optional promotion via a metadata_file.
    species_level = len(sciname.split()) >= 2
    if species_level and taxid:
        emit_taxid, emit_genus, emit_species = taxid, parts[0], parts[1]
    else:
        emit_taxid, emit_genus, emit_species = "", "NA", "NA"
        if sciname or taxid:
            note_bits.append(f"resolved_taxon={sciname or 'NA'}:{taxid or 'NA'}_not_species_level")

    return {
        "seq": seq,
        "upi": upi,
        "taxid": emit_taxid,
        "organism": sciname,
        "genus": emit_genus,
        "species": emit_species,
        "gene": gname or "NA",
        "description": pname or "NA",
        "n_organisms": n_orgs,
        "note": ";".join(note_bits),
    }


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
        out.append({
            "source_label": get(i_label) or default_label,
            "protein_id": get(i_pid) or fid,
            "fetch_id": fid,
            "source_accession": get(i_sacc) or acc or fid,
        })
    return out


def collect_requests(plan_path: Path) -> List[Dict[str, str]]:
    """Gather fetch requests from all uniparc mod_II rows."""
    out: List[Dict[str, str]] = []
    with plan_path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("enabled") or "").strip().lower() != "true":
                continue
            if (r.get("status") or "").strip().lower() == "disabled":
                continue
            if (r.get("module") or "").strip() != "mod_II":
                continue
            if parse_options(r.get("options") or "").get("fetch_source", "ncbi").lower() != "uniparc":
                continue
            label = (r.get("source_label") or "UNIPARC").strip()
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
    return _rows_to_requests(header, data, "UNIPARC")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", type=Path, help="Normalized source plan (selects uniparc mod_II rows)")
    ap.add_argument("-i", "--input", type=Path, help="Flat accession list (alternative to --plan)")
    ap.add_argument("-o", "--output", required=True, type=Path, help="FASTA out")
    ap.add_argument("--metadata", required=True, type=Path, help="user-metadata TSV out")
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--sleep", type=float, default=0.4)
    ap.add_argument("--retries", type=int, default=5)
    ap.add_argument("--batch-size", type=int, default=1, dest="batch_size",
                    help="Accepted for CLI symmetry; UniParc is fetched one UPI per request.")
    ap.add_argument("--organism-pick", choices=["best", "common", "first"], default="best",
                    dest="organism_pick",
                    help="Representative-organism rule for a UPI seen in multiple organisms.")
    ap.add_argument("--api-key", default="")
    ap.add_argument("--allow-missing", action="store_true",
                    help="Warn and skip accessions UniParc does not return instead of aborting.")
    ap.add_argument("--missing-out", type=Path, default=None,
                    help="Append skipped accessions to this TSV (source_label, accession, fetch_type, reason).")
    args = ap.parse_args()
    if bool(args.plan) == bool(args.input):
        die("provide exactly one of --plan or -i/--input")

    api_key = args.api_key or None
    requests = collect_requests(args.plan) if args.plan else collect_requests_from_list(args.input)
    if not requests:
        info("No UniParc accession rows in plan; nothing to do.")
        return

    seen = set()
    for r in requests:
        if r["protein_id"] in seen:
            die(f"duplicate protein_id across uniparc accession inputs: {r['protein_id']}")
        seen.add(r["protein_id"])

    info(f"UniParc accession requests: {len(requests)} (organism-pick={args.organism_pick})")

    resolved: Dict[str, dict] = {}   # fetch_id -> derived metadata dict
    for n, r in enumerate(requests, start=1):
        fid = r["fetch_id"]
        if not UPI_RE.match(fid):
            warn(f"'{fid}' does not look like a UniParc UPI; querying anyway.")
        entry = fetch_entry_json(fid, retries=args.retries, api_key=api_key)
        if entry is not None:
            meta = derive_metadata(entry, pick=args.organism_pick)
            if meta["seq"]:
                resolved[fid] = meta
            else:
                warn(f"{fid}: entry returned but no sequence; treating as missing.")
        if n % 50 == 0 or n == len(requests):
            info(f"Fetched {n}/{len(requests)}")
        time.sleep(args.sleep)

    missing_reqs = [r for r in requests if r["fetch_id"] not in resolved]
    kept_requests = [r for r in requests if r["fetch_id"] in resolved]
    missing_ids = [r["fetch_id"] for r in missing_reqs]

    if missing_ids:
        if not args.allow_missing:
            die(f"{len(missing_ids)} accession(s) not returned by UniParc: {missing_ids[:25]}"
                + (" ..." if len(missing_ids) > 25 else "")
                + "\n  Re-run with --allow-missing to skip these and continue.")
        warn(f"{len(missing_ids)} accession(s) not returned by UniParc; skipping (--allow-missing): "
             f"{missing_ids[:25]}" + (" ..." if len(missing_ids) > 25 else ""))
        if args.missing_out:
            args.missing_out.parent.mkdir(parents=True, exist_ok=True)
            new = not args.missing_out.exists()
            with args.missing_out.open("a", encoding="utf-8", newline="") as mo:
                if new:
                    mo.write("source_label\taccession\tfetch_type\treason\n")
                for r in missing_reqs:
                    mo.write(f"{r['source_label']}\t{r['fetch_id']}\tuniparc_accession\tnot_returned_by_uniparc\n")

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
            m = resolved[fid]
            seq = m["seq"]
            returned_id = m["upi"] or fid
            pname = m["description"] if m["description"] != "NA" else ""
            organism = m["organism"]

            # FASTA description: protein name + [organism] so build_db's stub
            # harvester can also recover the organism if metadata is ever lost.
            desc_bits = []
            if pname:
                desc_bits.append(pname)
            if organism:
                desc_bits.append(f"[{organism}]")
            fasta_desc = clean_desc(" ".join(desc_bits))

            out_header = (
                f">{clean_token(r['protein_id'])} "
                f"fetch_id={clean_token(fid)} "
                f"fetch_type=uniparc_accession "
                f"source_accession={clean_token(r['source_accession'])} "
                f"source_label={clean_token(r['source_label'])} "
                f"returned_fasta_id={clean_token(returned_id)}"
            )
            if fasta_desc:
                out_header += f" {fasta_desc}"
            fout.write(out_header + "\n")
            for chunk in wrap(seq):
                fout.write(chunk + "\n")

            mwr.writerow({
                "protein_id": clean_tsv(r["protein_id"]),
                "source": clean_tsv(r["source_label"]),
                # Real taxonomy derived from cross-references; blank taxid lets
                # build_db apply its 32644 default when UniParc had no organism.
                "ncbi_taxid": m["taxid"] or "",
                "genus": clean_tsv(m["genus"]),
                "species": clean_tsv(m["species"]),
                "strain_or_gene": clean_tsv(m["gene"]),
                "assembly": "NA",
                "description": clean_tsv(m["description"]),
                "fetch_id": clean_tsv(fid),
                "fetch_type": "uniparc_accession",
                "source_accession": clean_tsv(r["source_accession"]),
                "returned_fasta_id": clean_tsv(returned_id),
            })
            rwr.writerow({
                "source_label": r["source_label"], "protein_id": r["protein_id"],
                "fetch_id": fid, "fetch_type": "uniparc_accession",
                "source_accession": r["source_accession"], "returned_fasta_id": returned_id,
                "returned_description": (pname + (f" [{organism}]" if organism else "")).strip(),
                "sequence_length": len(seq),
                "status": "OK", "note": m["note"],
            })

        for r in missing_reqs:
            rwr.writerow({
                "source_label": r["source_label"], "protein_id": r["protein_id"],
                "fetch_id": r["fetch_id"], "fetch_type": "uniparc_accession",
                "source_accession": r["source_accession"], "returned_fasta_id": "NA",
                "returned_description": "", "sequence_length": 0,
                "status": "MISSING", "note": "not_returned_by_uniparc",
            })

    n = sum(1 for ln in args.output.open(encoding="utf-8") if ln.startswith(">"))
    if n != len(kept_requests):
        die(f"output record count mismatch: expected {len(kept_requests)}, wrote {n}")
    info(f"UniParc accession fetch complete: {n} sequences ({len(missing_ids)} missing)")
    info(f"FASTA: {args.output}")
    info(f"Metadata: {args.metadata}")
    info(f"Report: {args.report}")


if __name__ == "__main__":
    main()
