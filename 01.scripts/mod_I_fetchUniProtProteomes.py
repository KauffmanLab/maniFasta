#!/usr/bin/env python3

"""
mod_I_fetchUniProtProteomes.py
────────────────────────
maniFasta Module 1 (proteome / pan-proteome) — UniProt backend.

Parallels the NCBI path (mod_I_Fetch_UserSpecifiedSpecies_INFO.sh +
mod_I_download_UserSpecified_protein_FASTA.sh) but resolves and downloads whole
proteomes from UniProt instead of NCBI assemblies.

Two ways to select rows (mirrors mod_II's uniprot/uniparc/pdb fetchers):

  1. --plan   Whole-source declaration: reads the normalized source plan and
              processes mod_I rows where options contains fetch_source=uniprot.
              Per-row uniprot_proteome_type / uniprot_include_isoform /
              uniprot_reviewed options (if set) override the CLI defaults.

  2. -i/--input   Flat fetch list — an alternative to --plan for rows that
              came out of a *mixed* mod_I input_list, i.e. one whose own
              per-row `fetch_source` column selected uniprot for some rows
              and ncbi for others. mod_I_partition_proteome_sources.py
              splits such a list before this is called. Because these rows
              never had their own registry `options`, every row in the file
              shares one set of --proteome-type / --include-isoform /
              --reviewed flags.

Per input_list row, the UPID is resolved by priority:
    1. an explicit UPID  (column 'upid'/'proteome', or an Accession like UP000005640)
    2. taxid             -> proteomes search -> best proteome by type
    3. species name      -> proteomes search by organism_name -> best proteome

Output:
    <proteomes_dir>/<safe_species>__taxid<TAXID>__<UPID>.protein.faa  (one per row)
    appends rows to <proteomes_manifest> in the same schema the NCBI download
    writes, so build_db.py --proteomes_manifest consumes both transparently:

      source_label  Species_name  taxid  input_accession  picked_asm_or_acc
      source  status  out_faa  n_proteins  note

Stdlib only (urllib, csv, json, gzip).
"""

from __future__ import annotations

import argparse
import csv
import datetime
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

REST = "https://rest.uniprot.org"
TOOLNAME = "maniFasta_fetchUniProtProteomes"
UPID_RE = re.compile(r"^UP\d{9}$")

# ── Pan-proteome download (FTP) ───────────────────────────────────────────────
# The uniprotkb `proteome:UPID` stream returns ONLY that proteome's own set
# (e.g. the E. coli K-12 reference = ~4.4k proteins). The *pan* proteome — the
# non-redundant union across all strains of the species (~45k for E. coli) — is
# a pre-computed file published on the FTP site, NOT reachable by proteome:UPID.
#
# Pan proteomes are keyed on the *species* TAXID, not the UPID and not the
# taxonomic division. Layout (confirmed against the live FTP):
#   pan_proteomes/pp<taxid>/pp<taxid>.fasta.gz
# e.g. E. coli (taxid 562) -> pan_proteomes/pp562/pp562.fasta.gz
# See https://www.uniprot.org/help/pan_proteomes
PAN_FTP_BASE = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
                "knowledgebase/pan_proteomes")

# Higher rank = preferred when picking the "best" proteome for a taxon.
TYPE_RANK = {
    "reference proteome": 4,
    "representative proteome": 3,
    "other proteome": 2,
    "redundant proteome": 1,
}


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


def parse_bool(value: str, default: bool = False) -> bool:
    v = (value or "").strip().lower()
    if v in {"true", "t", "yes", "y", "1", "on"}:
        return True
    if v in {"false", "f", "no", "n", "0", "off", ""}:
        return False
    return default


def sanitize_name(s: str) -> str:
    s = re.sub(r"\s+", "_", (s or "").strip())
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def clean_one_line(s: str) -> str:
    return " ".join((s or "").replace("\r", " ").replace("\n", " ").split())


# ── HTTP ──────────────────────────────────────────────────────────────────────

def http_get(url: str, *, retries: int, sleep_sec: float,
             want_bytes: bool = False, api_key: Optional[str] = None):
    """GET with retries + Retry-After handling. Returns text (or bytes)."""
    headers = {"User-Agent": TOOLNAME}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"  # harmless if unsupported
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=300) as resp:
                raw = resp.read()
                return raw if want_bytes else raw.decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 429:
                wait = float(e.headers.get("Retry-After", 2 ** attempt))
                warn(f"429 rate limited; sleeping {wait}s")
                time.sleep(wait)
                continue
            warn(f"HTTP {e.code} (attempt {attempt}/{retries}): {url}")
        except Exception as e:  # noqa: BLE001
            last_err = e
            warn(f"GET attempt {attempt}/{retries} failed: {e}")
        if attempt < retries:
            time.sleep(min(2 ** attempt, 30))
    raise RuntimeError(f"GET failed after {retries} attempts ({url}): {last_err}")


def http_get_optional(url: str, *, retries: int, want_bytes: bool = False,
                      api_key: Optional[str] = None):
    """GET that treats 404/410 as a clean miss (returns None) instead of raising.

    Used when probing candidate FTP pan-proteome URLs: a missing candidate must
    fail fast (no exponential-backoff retry storm on an expected 404), while a
    real network hiccup (timeouts, 5xx, 429) still gets a couple of retries.
    """
    headers = {"User-Agent": TOOLNAME}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=300) as resp:
                raw = resp.read()
                return raw if want_bytes else raw.decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code in (404, 410):
                return None  # expected miss — do not retry
            last_err = e
            if e.code == 429:
                wait = float(e.headers.get("Retry-After", 2 ** attempt))
                warn(f"429 rate limited; sleeping {wait}s")
                time.sleep(wait)
                continue
            warn(f"HTTP {e.code} (attempt {attempt}/{retries}): {url}")
        except Exception as e:  # noqa: BLE001
            last_err = e
            warn(f"GET attempt {attempt}/{retries} failed: {e}")
        if attempt < retries:
            time.sleep(min(2 ** attempt, 30))
    if last_err is not None:
        warn(f"pan probe GET gave up on {url}: {last_err}")
    return None


# ── UniProt lookups ─────────────────────────────────────────────────────────

def proteomes_search(query: str, *, retries: int, sleep_sec: float,
                     api_key: Optional[str]) -> List[dict]:
    params = {"query": query, "format": "json", "size": "50"}
    url = f"{REST}/proteomes/search?{urllib.parse.urlencode(params)}"
    txt = http_get(url, retries=retries, sleep_sec=sleep_sec, api_key=api_key)
    try:
        return json.loads(txt).get("results", []) or []
    except json.JSONDecodeError as e:
        warn(f"proteomes search JSON parse error: {e}")
        return []


def proteome_entry(upid: str, *, retries: int, sleep_sec: float,
                   api_key: Optional[str]) -> dict:
    url = f"{REST}/proteomes/{upid}?format=json"
    txt = http_get(url, retries=retries, sleep_sec=sleep_sec, api_key=api_key)
    try:
        return json.loads(txt)
    except json.JSONDecodeError:
        return {}


def _result_fields(r: dict) -> tuple[str, str, str, str]:
    """Return (upid, proteome_type, taxid, organism) from a proteomes result."""
    upid = (r.get("id") or r.get("upid") or "").strip()
    ptype = (r.get("proteomeType") or "").strip()
    tax = r.get("taxonomy") or {}
    taxid = str(tax.get("taxonId") or "").strip()
    organism = (tax.get("scientificName") or "").strip()
    return upid, ptype, taxid, organism


def pick_best(results: List[dict], want: str) -> Optional[dict]:
    if not results:
        return None
    # Drop Redundant/Excluded proteomes: these records still show up in
    # proteomes_search but their sequences aren't downloadable (a side effect
    # of the 2015 UniProt proteome redundancy minimization). Without this
    # filter pick_best can hand back a stub UPID that resolves cleanly and
    # then produces an empty FASTA at download time -- exactly the "empty
    # proteome for UP000030883" case seen on Lactobacillus acidophilus.
    # If EVERY result is a stub, fall back to the original unfiltered list
    # so we don't turn a bad-pick bug into a spurious no_proteome_found.
    good = [r for r in results
            if (r.get("proteomeType") or "").strip().lower()
            not in ("redundant proteome", "excluded proteome")]
    pool = good or results
    if want in ("reference", "representative"):
        target = f"{want} proteome"
        exact = [r for r in pool
                 if (r.get("proteomeType") or "").strip().lower() == target]
        if exact:
            return exact[0]
        # fall through to best-by-rank if the exact type is unavailable
    if want == "any":
        return pool[0]
    return max(pool, key=lambda r: TYPE_RANK.get(
        (r.get("proteomeType") or "").strip().lower(), 0))


def pan_upid_from_entry(entry: dict) -> Optional[str]:
    """Return the linked pan-proteome UPID from an already-fetched proteome entry."""
    for key in ("panProteome", "panproteome", "panProteomeUpid"):
        val = entry.get(key)
        if isinstance(val, str) and UPID_RE.match(val.strip()):
            return val.strip()
        if isinstance(val, dict):
            cand = (val.get("upid") or val.get("id") or "").strip()
            if UPID_RE.match(cand):
                return cand
    return None


def resolve_pan_upid(upid: str, *, retries: int, sleep_sec: float,
                     api_key: Optional[str]) -> Optional[str]:
    """Read the proteome entry and return its linked pan-proteome UPID, if any."""
    entry = proteome_entry(upid, retries=retries, sleep_sec=sleep_sec, api_key=api_key)
    return pan_upid_from_entry(entry)


# ── Pan-proteome (FTP union) ─────────────────────────────────────────────────

_HREF_RE = re.compile(r'href=["\']?([^"\'>\s]+\.fasta\.gz)', re.IGNORECASE)


def _decode_fasta_bytes(raw: bytes) -> str:
    try:
        return gzip.decompress(raw).decode("utf-8", "replace")
    except (OSError, EOFError):
        return raw.decode("utf-8", "replace")  # server delivered plain text


def species_taxid_from_entry(entry: Optional[dict]) -> str:
    """Return the species-rank taxId from a proteome entry lineage, if present.

    Lets a strain-level input taxid still resolve to the species-keyed pan file
    (pan proteomes are grouped at the species level).
    """
    if not isinstance(entry, dict):
        return ""
    lineage = entry.get("taxonLineage")
    if isinstance(lineage, list):
        for node in lineage:
            if isinstance(node, dict) and str(node.get("rank", "")).strip().lower() == "species":
                tid = str(node.get("taxonId") or "").strip()
                if tid.isdigit():
                    return tid
    return ""


def _pan_taxid_candidates(taxid: str, entry: Optional[dict]) -> List[str]:
    """Ordered, de-duplicated numeric taxids to try as pp<taxid> (row taxid first,
    then the species-rank taxid from the entry for strain-level inputs)."""
    out: List[str] = []
    for t in ((taxid or "").strip(), species_taxid_from_entry(entry)):
        if t.isdigit() and t not in out:
            out.append(t)
    return out


def download_pan_proteome_fasta(base_upid: str, pan_upid: Optional[str],
                                taxid: str, entry: Optional[dict], *,
                                retries: int, sleep_sec: float,
                                api_key: Optional[str]) -> tuple[str, str, str]:
    """Fetch the pre-computed pan-proteome FASTA (non-redundant union across the
    species' strains) from the UniProt FTP.

    Pan proteomes are keyed on the *species* taxid:
        pan_proteomes/pp<taxid>/pp<taxid>.fasta.gz

    Returns (fasta_text, source_note, pan_id). fasta_text is "" if no pan file
    exists for this taxon (caller then falls back to the single reference
    proteome). pan_id is e.g. "pp562" on success, "" otherwise.
    """
    cands = _pan_taxid_candidates(taxid, entry)
    if not cands:
        return "", "pan_no_taxid", ""

    for t in cands:
        folder = f"{PAN_FTP_BASE}/pp{t}"
        url = f"{folder}/pp{t}.fasta.gz"
        raw = http_get_optional(url, retries=retries, want_bytes=True, api_key=api_key)
        if raw:
            time.sleep(sleep_sec)
            return _decode_fasta_bytes(raw), f"pan_ftp:{url}", f"pp{t}"

        # Fallback: folder exists but the protein FASTA is named differently.
        # List pp<taxid>/ and prefer pp<taxid>.fasta.gz, else the first plain
        # protein *.fasta.gz (skip DNA/additional sidecars).
        idx = http_get_optional(f"{folder}/", retries=retries, api_key=api_key)
        if idx:
            names = [h.rsplit("/", 1)[-1] for h in _HREF_RE.findall(idx)]
            pick = next((n for n in names if n == f"pp{t}.fasta.gz"), None)
            if not pick:
                pick = next((n for n in names
                             if n.endswith(".fasta.gz")
                             and "_DNA" not in n and "additional" not in n.lower()), None)
            if pick:
                raw = http_get_optional(f"{folder}/{pick}",
                                        retries=retries, want_bytes=True, api_key=api_key)
                if raw:
                    time.sleep(sleep_sec)
                    return _decode_fasta_bytes(raw), f"pan_ftp:{folder}/{pick}", f"pp{t}"

    return "", "pan_not_published_for_taxon", ""


def resolve_row(species: str, taxid: str, accession: str, upid_col: str,
                want_type: str, *, retries: int, sleep_sec: float,
                api_key: Optional[str]) -> dict:
    """Resolve one input row to {upid, taxid, organism, note}. upid empty on failure."""
    explicit = (upid_col or "").strip()
    if not explicit and UPID_RE.match((accession or "").strip()):
        explicit = accession.strip()

    chosen = None
    note_bits: List[str] = []

    base_entry: Optional[dict] = None
    if explicit and UPID_RE.match(explicit):
        entry = proteome_entry(explicit, retries=retries, sleep_sec=sleep_sec, api_key=api_key)
        time.sleep(sleep_sec)
        base_entry = entry
        _, ptype, etax, eorg = _result_fields(entry)
        chosen = {"upid": explicit, "proteomeType": ptype,
                  "taxonomy": {"taxonId": etax or taxid, "scientificName": eorg or species}}
        note_bits.append("explicit_upid")
    elif (taxid or "").strip():
        # organism_id matches the EXACT node the proteome is filed under
        # (species or strain). This is right for species-level taxIDs and
        # also for keeping the current single-species behavior when it works.
        # If it comes back empty, retry with taxonomy_id, which is
        # lineage-aware: it matches any proteome whose organism sits at or
        # under the given node. That's what lets higher-rank inputs
        # (family, genus, order) resolve at all -- previously they died
        # with no_proteome_found because organism_id=<family> matches no
        # proteome record directly. pick_best still collapses the lineage
        # set to a single winner, so per-row output count is unchanged
        # (one FASTA per row); clade expansion is a separate future feature.
        t = taxid.strip()
        results = proteomes_search(f"(organism_id:{t})",
                                   retries=retries, sleep_sec=sleep_sec, api_key=api_key)
        time.sleep(sleep_sec)
        scope = "by_taxid"
        if not results:
            results = proteomes_search(f"(taxonomy_id:{t})",
                                       retries=retries, sleep_sec=sleep_sec, api_key=api_key)
            time.sleep(sleep_sec)
            scope = "by_taxid_lineage"
        chosen = pick_best(results, want_type)
        note_bits.append(f"{scope}:{len(results)}_candidates")
    elif (species or "").strip():
        q = f'(organism_name:"{species.strip()}")'
        results = proteomes_search(q, retries=retries, sleep_sec=sleep_sec, api_key=api_key)
        time.sleep(sleep_sec)
        chosen = pick_best(results, want_type)
        note_bits.append(f"by_name:{len(results)}_candidates")
    else:
        return {"upid": "", "taxid": "", "organism": "", "note": "no_identifier"}

    if not chosen:
        return {"upid": "", "taxid": taxid, "organism": species, "note": "no_proteome_found"}

    upid, ptype, rtax, rorg = _result_fields(chosen)
    if not upid:
        return {"upid": "", "taxid": taxid, "organism": species, "note": "result_missing_upid"}

    pan_upid: Optional[str] = None
    pan_entry: Optional[dict] = None
    if want_type == "pan":
        # Reuse the entry we already fetched for an explicit UPID; otherwise
        # fetch the chosen proteome's entry once (needed for the panProteome
        # link and for taxonomic-division detection on the FTP).
        if base_entry is not None and (base_entry.get("id") or "").strip() == upid:
            pan_entry = base_entry
        else:
            pan_entry = proteome_entry(upid, retries=retries, sleep_sec=sleep_sec, api_key=api_key)
            time.sleep(sleep_sec)
        pan_upid = pan_upid_from_entry(pan_entry or {})
        note_bits.append(f"pan_of:{upid}" if pan_upid else "no_pan_link_probe_ftp_on_ref")

    return {
        "upid": upid,                 # base reference/best proteome UPID (unchanged)
        "pan_upid": pan_upid or "",   # linked pan UPID if UniProt exposes one
        "pan_entry": pan_entry,       # proteome entry (for FTP division detection)
        "taxid": (rtax or taxid or "").strip(),
        "organism": (rorg or species or "").strip(),
        "note": ";".join([p for p in [ptype, *note_bits] if p]),
    }


def download_proteome_fasta(upid: str, *, include_isoform: bool, reviewed: bool,
                            retries: int, sleep_sec: float,
                            api_key: Optional[str]) -> str:
    query = f"(proteome:{upid})"
    if reviewed:
        query = f"{query} AND (reviewed:true)"
    params = {
        "query": query,
        "format": "fasta",
        "compressed": "true",
        "includeIsoform": "true" if include_isoform else "false",
    }
    url = f"{REST}/uniprotkb/stream?{urllib.parse.urlencode(params)}"
    raw = http_get(url, retries=retries, sleep_sec=sleep_sec, want_bytes=True, api_key=api_key)
    try:
        return gzip.decompress(raw).decode("utf-8", "replace")
    except (OSError, EOFError):
        return raw.decode("utf-8", "replace")  # server ignored compression


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


def iter_uniprot_proteome_rows(plan_path: Path):
    """Yield (source_label, options, input_list_path) for uniprot proteome rows."""
    with plan_path.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if (r.get("enabled") or "").strip().lower() != "true":
                continue
            if (r.get("status") or "").strip().lower() == "disabled":
                continue
            if (r.get("module") or "").strip() != "mod_I":
                continue
            opts = parse_options(r.get("options") or "")
            if opts.get("fetch_source", "ncbi").lower() != "uniprot":
                continue
            yield (r.get("source_label") or "UNIPROT").strip(), opts, (r.get("input_list") or "").strip()


def read_flat_list_rows(list_path: Path) -> List[dict]:
    """Read a flat, already-partitioned proteome fetch list (as written by
    mod_I_partition_proteome_sources.py from the per-row `fetch_source`
    column) into row dicts. Mirrors the mod_II *_fromlist readers: no
    per-row `options` here, since options live on the registry row, not on
    an individual accession/species line -- callers use the -i mode's CLI
    flags (--proteome-type / --include-isoform / --reviewed) as the shared
    default for every row in the file.
    """
    header, data = read_clean_table(list_path)
    i_label = col_index(header, "source_label")
    i_spp = col_index(header, "species", "Species name", "species name")
    i_tax = col_index(header, "taxonID", "TaxonID", "taxid")
    i_acc = col_index(header, "accession", "Accession")
    i_up = col_index(header, "upid", "proteome", "UPID")
    if i_spp == -1 and i_acc == -1 and i_tax == -1:
        die(f"{list_path}: needs a 'species', 'accession', or 'taxonID' column. Found: {header}")

    def get(row: List[str], idx: int) -> str:
        return (row[idx].strip() if 0 <= idx < len(row) else "")

    rows = []
    for row in data:
        species = clean_one_line(get(row, i_spp))
        taxid_in = get(row, i_tax).replace(" ", "")
        accession = get(row, i_acc)
        upid_col = get(row, i_up)
        if not (species or taxid_in or accession or upid_col):
            continue
        rows.append({
            "source_label": get(row, i_label) or "UNIPROT",
            "species": species,
            "taxid_in": taxid_in,
            "accession": accession,
            "upid_col": upid_col,
        })
    return rows


# ── Main ────────────────────────────────────────────────────────────────────

def fetch_one_row(*, label: str, species: str, taxid_in: str, accession: str,
                  upid_col: str, want_type: str, include_isoform: bool,
                  reviewed: bool, args, api_key: Optional[str],
                  mf, counts: Dict[str, int]) -> None:
    """Resolve one row to a UniProt proteome, download it, and write one
    manifest line. Shared by both --plan mode (per-registry-row options) and
    -i/--input flat-list mode (shared CLI-level options) so the two entry
    points can never drift apart."""
    src_tag = "uniprot_pan" if want_type == "pan" else "uniprot"
    disp = species or accession or upid_col or taxid_in or "<row>"

    res = resolve_row(species, taxid_in, accession, upid_col, want_type,
                      retries=args.retries, sleep_sec=args.sleep, api_key=api_key)
    upid = res["upid"]
    taxid = res["taxid"] or taxid_in
    organism = res["organism"] or species or "NA"

    # A pan request only needs the species taxid to locate pp<taxid> on the FTP,
    # so don't reject it just because UPID resolution came up empty (e.g. a
    # taxid-only row). Non-pan rows still require a resolved proteome UPID.
    if not upid and not (want_type == "pan" and (taxid or "").strip().isdigit()):
        warn(f"  [{disp}] resolve failed: {res['note']}")
        mf.write("\t".join([label, organism, taxid, accession or upid_col, "",
                            src_tag, "FAILED_RESOLVE", "", "0", res["note"]]) + "\n")
        counts["fail"] += 1
        return

    picked = upid          # what lands in the manifest 'picked_asm_or_acc' column
    pan_source_note = ""
    fasta: Optional[str] = None

    if want_type == "pan":
        # Pan proteome = the non-redundant union across all strains of the
        # species (~45k for E. coli). That set is NOT reachable via a
        # `proteome:UPID` uniprotkb query -- that only returns the single
        # reference proteome (~4.4k). The pan set is a pre-computed FTP file
        # keyed on the species taxid (pan_proteomes/pp<taxid>/); fetch that.
        pan_upid = (res.get("pan_upid") or "").strip()
        pan_entry = res.get("pan_entry")
        pan_id = ""
        try:
            fasta, pan_source_note, pan_id = download_pan_proteome_fasta(
                upid, pan_upid or None, taxid, pan_entry,
                retries=args.retries, sleep_sec=args.sleep, api_key=api_key)
        except RuntimeError as e:
            warn(f"  [{disp}] pan download error for taxid {taxid or '?'}: {e}")
            fasta, pan_source_note, pan_id = "", f"pan_error:{str(e)[:80]}", ""

        if fasta:
            picked = pan_id or (pan_upid or upid)   # e.g. "pp562"
        else:
            # No pan file published for this taxon (or the FTP was unreachable).
            # Fall back to the single reference proteome if we have one; if this
            # was a taxid-only row with no resolved UPID, there is nothing to
            # fall back to, so fail the row clearly.
            if not upid:
                warn(f"  [{disp}] no pan file and no reference proteome resolved "
                     f"({pan_source_note})")
                mf.write("\t".join([label, organism, taxid, accession or upid_col, "",
                                    src_tag, "FAILED_RESOLVE", "", "0",
                                    f"{res['note']};{pan_source_note}"]) + "\n")
                counts["fail"] += 1
                time.sleep(args.sleep)
                return
            warn(f"  [{disp}] no pan-proteome file found ({pan_source_note}); "
                 f"falling back to reference proteome {upid}")
            try:
                fasta = download_proteome_fasta(
                    upid, include_isoform=include_isoform, reviewed=reviewed,
                    retries=args.retries, sleep_sec=args.sleep, api_key=api_key)
            except RuntimeError as e:
                warn(f"  [{disp}] fallback download failed for {upid}: {e}")
                mf.write("\t".join([label, organism, taxid, accession or upid_col, upid,
                                    src_tag, "FAILED_DOWNLOAD", "", "0", str(e)[:120]]) + "\n")
                counts["fail"] += 1
                time.sleep(args.sleep)
                return
    else:
        try:
            fasta = download_proteome_fasta(
                upid, include_isoform=include_isoform, reviewed=reviewed,
                retries=args.retries, sleep_sec=args.sleep, api_key=api_key)
        except RuntimeError as e:
            warn(f"  [{disp}] download failed for {upid}: {e}")
            mf.write("\t".join([label, organism, taxid, accession or upid_col, upid,
                                src_tag, "FAILED_DOWNLOAD", "", "0", str(e)[:120]]) + "\n")
            counts["fail"] += 1
            time.sleep(args.sleep)
            return

    n_prot = fasta.count("\n>") + (1 if fasta.startswith(">") else 0)
    if n_prot == 0 or not fasta.lstrip().startswith(">"):
        warn(f"  [{disp}] empty proteome for {picked}")
        mf.write("\t".join([label, organism, taxid, accession or upid_col, picked,
                            src_tag, "FAILED_EMPTY", "", "0", "no sequences returned"]) + "\n")
        counts["fail"] += 1
        time.sleep(args.sleep)
        return

    # Compose the final note (pan path annotates the source it used). Status
    # stays exactly "OK" so build_db.py's `status != "OK"` filter keeps the row;
    # the fallback is made unmistakable via the source column (uniprot_pan) and
    # the note (grep 'pan_unavailable' to audit which taxa lacked a pan file).
    note = res["note"]
    status = "OK"
    if want_type == "pan":
        if pan_source_note.startswith("pan_ftp:"):
            note = f"{note};{pan_source_note}"
        else:
            note = f"{note};pan_unavailable_fell_back_to_reference:{pan_source_note}"

    safe = sanitize_name(organism if organism != "NA" else picked)
    out = args.proteomes_dir / f"{safe}__taxid{taxid or 'NA'}__{picked}.protein.faa"
    out.write_text(fasta if fasta.endswith("\n") else fasta + "\n", encoding="utf-8")

    info(f"  [{disp}] -> {picked} ({note}) proteins={n_prot}")
    mf.write("\t".join([label, organism, taxid, accession or upid_col, picked,
                        src_tag, status, str(out), str(n_prot), note]) + "\n")
    mf.flush()
    counts["ok"] += 1
    time.sleep(args.sleep)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", type=Path,
                    help="Normalized source plan (selects whole-source mod_I "
                         "rows with options fetch_source=uniprot)")
    ap.add_argument("-i", "--input", type=Path,
                    help="Flat proteome fetch list -- alternative to --plan. "
                         "Rows come from a mixed mod_I input_list whose "
                         "per-row fetch_source column selected uniprot "
                         "(as partitioned by mod_I_partition_proteome_sources.py). "
                         "Every row in the file shares --proteome-type / "
                         "--include-isoform / --reviewed since options here "
                         "are CLI-level, not per-registry-row.")
    ap.add_argument("--proteomes-dir", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path,
                    help="proteomes_manifest.tsv (created or appended)")
    ap.add_argument("--proteome-type", default="best", dest="proteome_type",
                    choices=["best", "reference", "representative", "any", "pan"],
                    help="Global default proteome type when a source row does not set "
                         "uniprot_proteome_type in its options (default best). A per-source "
                         "uniprot_proteome_type option still overrides this in --plan mode.")
    ap.add_argument("--include-isoform", action="store_true", dest="include_isoform",
                    help="(-i mode only) include isoforms; --plan mode reads this "
                         "per-row from uniprot_include_isoform instead.")
    ap.add_argument("--reviewed", action="store_true",
                    help="(-i mode only) Swiss-Prot only; --plan mode reads this "
                         "per-row from uniprot_reviewed instead.")
    ap.add_argument("--sleep", type=float, default=0.4,
                    help="seconds between UniProt requests (default 0.4)")
    ap.add_argument("--retries", type=int, default=5)
    ap.add_argument("--api-key", default="", help="optional UniProt token")
    args = ap.parse_args()

    if bool(args.plan) == bool(args.input):
        die("provide exactly one of --plan or -i/--input")

    api_key = args.api_key or None
    args.proteomes_dir.mkdir(parents=True, exist_ok=True)

    manifest_cols = ["source_label", "Species_name", "taxid", "input_accession",
                     "picked_asm_or_acc", "source", "status", "out_faa",
                     "n_proteins", "note"]
    new_manifest = not args.manifest.exists()
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    mf = args.manifest.open("a", encoding="utf-8", newline="")
    if new_manifest:
        mf.write("\t".join(manifest_cols) + "\n")

    counts = {"ok": 0, "fail": 0}

    if args.plan:
        rows = list(iter_uniprot_proteome_rows(args.plan))
        if not rows:
            info("No UniProt proteome rows in plan; nothing to do.")
            mf.close()
            return

        for label, opts, input_list in rows:
            want_type = opts.get("uniprot_proteome_type", args.proteome_type).lower()
            include_isoform = parse_bool(opts.get("uniprot_include_isoform", ""), default=False)
            reviewed = parse_bool(opts.get("uniprot_reviewed", ""), default=False)

            src = Path(input_list)
            if not src.is_file():
                die(f"input_list not found for '{label}': {src}")

            header, data = read_clean_table(src)
            i_spp = col_index(header, "species", "Species name", "species name")
            i_tax = col_index(header, "taxonID", "TaxonID", "taxid")
            i_acc = col_index(header, "accession", "Accession")
            i_up = col_index(header, "upid", "proteome", "UPID")

            info(f"[{label}] {len(data)} row(s); type={want_type} "
                 f"isoform={include_isoform} reviewed={reviewed}")

            for row in data:
                def get(idx: int) -> str:
                    return (row[idx].strip() if 0 <= idx < len(row) else "")

                fetch_one_row(
                    label=label,
                    species=clean_one_line(get(i_spp)),
                    taxid_in=get(i_tax).replace(" ", ""),
                    accession=get(i_acc),
                    upid_col=get(i_up),
                    want_type=want_type, include_isoform=include_isoform,
                    reviewed=reviewed, args=args, api_key=api_key,
                    mf=mf, counts=counts,
                )
    else:
        flat_rows = read_flat_list_rows(args.input)
        if not flat_rows:
            info("No UniProt proteome rows in input list; nothing to do.")
            mf.close()
            return

        want_type = args.proteome_type.lower()
        info(f"[flat-list] {len(flat_rows)} row(s); type={want_type} "
             f"isoform={args.include_isoform} reviewed={args.reviewed}")

        for r in flat_rows:
            fetch_one_row(
                label=r["source_label"], species=r["species"],
                taxid_in=r["taxid_in"], accession=r["accession"],
                upid_col=r["upid_col"], want_type=want_type,
                include_isoform=args.include_isoform, reviewed=args.reviewed,
                args=args, api_key=api_key, mf=mf, counts=counts,
            )

    mf.close()
    info(f"UniProt proteome fetch complete: OK={counts['ok']}, failed={counts['fail']}")
    info(f"Manifest: {args.manifest}")
    if counts["ok"] == 0:
        die("no UniProt proteomes were downloaded successfully")


if __name__ == "__main__":
    main()
