#!/usr/bin/env python3

import argparse
from dataclasses import dataclass, field, fields
from pathlib import Path
import re
import csv
import sys
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

# Python csv has a conservative default field-size limit (~128 KB).
# Some source FASTA descriptions / metadata notes can be longer than that,
# especially after preserving full provenance strings. Raise the limit once
# at import time so DictReader can parse long description fields safely.
_csv_limit = sys.maxsize
while True:
    try:
        csv.field_size_limit(_csv_limit)
        break
    except OverflowError:
        _csv_limit = int(_csv_limit / 10)

DB_PREFIX    = "maniFastaDB"
EUTILS_BASE  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EUTILS_TOOL  = "maniFastaDB_build"

# -----------------------------------------------------------
# FASTA ITERATION UTILITIES
# -----------------------------------------------------------

def iter_fasta(path: Path):
    """Yield (header, seq) tuples from a FASTA file."""
    header = None
    seq_chunks = []

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())

        if header is not None:
            yield header, "".join(seq_chunks)


def wrap_sequence(seq: str, width: int = 60):
    """Wrap sequence lines to width."""
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


# -----------------------------------------------------------
# STRICT VALIDATORS
# -----------------------------------------------------------

def require_taxid(taxid: str, context: str):
    taxid = (taxid or "").strip()

    if not taxid or taxid == "NA":
        raise RuntimeError(f"HARD FAIL: missing taxid for {context}")

    if not re.fullmatch(r"[0-9]+", taxid):
        raise RuntimeError(f"HARD FAIL: non-numeric taxid '{taxid}' for {context}")

    if taxid.lstrip("0") == "":
        raise RuntimeError(
            f"HARD FAIL: placeholder taxid '{taxid}' for {context}. "
            f"eHOMD writes 0 in 'NCBI Taxon ID' for taxa it has not yet mapped. "
            f"Re-fetch the eHOMD taxon table, or supply an override."
        )

    return taxid

# eHOMD changed HMT-ID zero-padding between releases (HMT-047 -> HMT-0047,
# eHOMD changed HMT-ID zero-padding between releases (HMT-047 -> HMT-0047,
# observed 2026-07-31 vs 2026-08-02) while PROKKA/V11.03/GCA_ID_info.txt kept
# the old 3-digit form. Canonicalise to 4-digit everywhere -- join keys and
# emitted HMT= / hmt_id values -- so output is stable across upstream changes.
_HMT_RE = re.compile(r"^HMT-?0*([0-9]+)$", re.IGNORECASE)

def norm_hmt(hmt: str) -> str:
    hmt = (hmt or "").strip()
    m = _HMT_RE.match(hmt)
    return f"HMT-{int(m.group(1)):04d}" if m else hmt

def require_organism(org: str, context: str):
    org = (org or "").strip()

    if not org or org == "NA":
        raise RuntimeError(f"HARD FAIL: missing organism for {context}")

    return org

def safe_one_line(s: str):
    return re.sub(r"[\t\r\n]+", " ", (s or "")).strip()


def _ts() -> str:
    import datetime
    return datetime.datetime.now().strftime("%H:%M:%S")


def info(msg: str):
    print(f"[INFO] [{_ts()}] {msg}", file=sys.stderr)


def warn(msg: str):
    print(f"[WARN] [{_ts()}] {msg}", file=sys.stderr)


# -----------------------------------------------------------
# LOAD METADATA TABLES: HOMD
# -----------------------------------------------------------

def load_gca_to_hmt(path: Path):
    """
    Parse GCA_ID_info.txt.

    Expected whitespace-separated columns:
      GCA-ID  HMT-ID  Genus  Species  Strain  Contigs  Combined_Size  Sequence_Source
    """
    mapping = {}

    with path.open() as fh:
        for line in fh:
            line = line.strip()

            if not line or line.startswith("GCA-ID"):
                continue

            cols = line.split()

            if len(cols) < 5:
                continue

            gca_id  = cols[0]
            hmt_id  = cols[1]
            genus   = cols[2]
            species = cols[3]

            if len(cols) > 8:
                strain_tokens = cols[4:-3]
                strain = " ".join(strain_tokens)
            else:
                strain = cols[4]

            mapping[gca_id] = {
                "hmt":     hmt_id,
                "genus":   genus,
                "species": species,
                "strain":  strain if strain else "NA",
            }

    if not mapping:
        raise RuntimeError(f"HARD FAIL: No entries parsed from {path}; check delimiter/format.")

    return mapping


def load_hmt_to_tax(path: Path):
    """
    Parse HOMD taxonomy table (full or body-site-filtered by homd_filter.py).

    Expects tab-delimited table with a header containing:
      HMT-ID, Genus, Species, NCBI Taxon ID

    When a body-site-filtered TSV is supplied the HMT IDs present in this
    table define which proteins are included; proteins whose HMT is absent
    are silently filtered out during the FASTA processing step.
    """
    mapping = {}
    header  = None
    idx     = {}

    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("HOMD.org Taxon Data::"):
                continue

            if header is None and line.startswith("HMT-ID"):
                header = line.split("\t")
                idx    = {name: i for i, name in enumerate(header)}

                need = {"HMT-ID", "NCBI Taxon ID", "Genus", "Species"}

                if not need.issubset(idx.keys()):
                    raise RuntimeError(f"HARD FAIL: Unexpected HOMD tax header: {header}")

                continue

            if header is None:
                continue

            cols = line.split("\t")

            if len(cols) < len(header):
                cols += [""] * (len(header) - len(cols))

            hmt_id = cols[idx["HMT-ID"]].strip()

            if not hmt_id:
                continue

            taxid   = cols[idx["NCBI Taxon ID"]].strip() or "NA"
            genus   = cols[idx["Genus"]].strip()          or "NA"
            species = cols[idx["Species"]].strip()        or "NA"

            mapping[norm_hmt(hmt_id)] = {
                "taxid":   taxid,
                "genus":   genus,
                "species": species,
            }

    if not mapping:
        raise RuntimeError(
            f"HARD FAIL: No HMT-to-taxid entries parsed from {path}. "
            f"If a body-site filter was applied, check that the requested "
            f"site(s) exist in the HOMD taxonomy file."
        )

    return mapping


# -----------------------------------------------------------
# UNIPROT-LIKE HEADER PARSER
# HUMAN / SARS2 / CRAP
# -----------------------------------------------------------

def parse_uniprot_like_header(header: str):
    """
    Parse UniProt-style headers that contain OS= and OX=.
    HARD FAIL if OX or OS missing.
    """
    parts   = header.split()
    id_part = parts[0] if parts else ""
    desc    = " ".join(parts[1:]) if len(parts) > 1 else ""

    tokens  = id_part.split("|")
    uniprot = tokens[1] if len(tokens) > 1 else tokens[0] if tokens else "NA"

    m_gene    = re.search(r"\bGN=([A-Za-z0-9_.-]+)", header)
    gene      = m_gene.group(1) if m_gene else "NA"

    m_tax     = re.search(r"\bOX=([0-9]+)", header)
    taxid     = m_tax.group(1) if m_tax else "NA"

    m_species = re.search(r"\bOS=([^=]+?)(?:\s[A-Z]{2}=|$)", header)
    organism  = m_species.group(1).strip() if m_species else "NA"

    taxid    = require_taxid(taxid, f"uniprot_like:{uniprot}")
    organism = require_organism(organism, f"uniprot_like:{uniprot}")

    return {
        "id":          uniprot,
        "gene":        gene,
        "taxid":       taxid,
        "organism":    organism,
        "description": desc,
    }


def split_organism_for_map(organism: str, fallback_genus: str = "NA"):
    toks = organism.split()

    if not toks:
        return fallback_genus, "NA"

    genus   = toks[0]
    species = " ".join(toks[1:]) if len(toks) > 1 else organism

    return genus, species


# -----------------------------------------------------------
# PROTEOMES MANIFEST LOADER
# -----------------------------------------------------------

def load_proteomes_manifest(manifest_path: Path):
    """
    Reads proteomes_manifest.tsv.

    Supports both older and newer manifest schemas.

    Old expected columns:
      Species_name, taxid, picked_assembly, assembly_source, status, out_faa

    New expected columns:
      Species_name, taxid, input_accession, picked_asm_or_acc, source,
      status, out_faa, n_proteins, note

    Internally always returns:
      Species_name, taxid, picked_assembly, assembly_source, out_faa
    """
    rows                 = []
    skipped_non_ok       = 0
    skipped_missing_file = 0
    failed_rows: list    = []   # [{species, reason}] surfaced in build summary

    with manifest_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        if reader.fieldnames is None:
            raise RuntimeError(f"HARD FAIL: {manifest_path} has no header row.")

        fieldnames = set(reader.fieldnames)

        old_schema = {
            "Species_name", "taxid", "picked_assembly",
            "assembly_source", "status", "out_faa",
        }
        new_schema = {
            "Species_name", "taxid", "picked_asm_or_acc",
            "source", "status", "out_faa",
        }
        # source_label is optional (absent in older manifests)
        has_source_label = "source_label" in fieldnames

        if old_schema.issubset(fieldnames):
            schema     = "old"
            picked_col = "picked_assembly"
            source_col = "assembly_source"
        elif new_schema.issubset(fieldnames):
            schema     = "new"
            picked_col = "picked_asm_or_acc"
            source_col = "source"
        else:
            raise RuntimeError(
                "HARD FAIL: manifest does not match expected old or new schema.\n"
                f"Found columns: {reader.fieldnames}\n"
                f"Old required columns: {sorted(old_schema)}\n"
                f"New required columns: {sorted(new_schema)}"
            )

        info(f"Manifest schema detected: {schema}")

        for r in reader:
            spp    = (r.get("Species_name") or "").strip()
            status = (r.get("status")       or "").strip()

            if status != "OK":
                skipped_non_ok += 1
                warn(f"Skipping manifest row because status={status}: {spp}")
                failed_rows.append({"species": spp, "reason": f"status={status}"})
                continue

            taxid = require_taxid((r.get("taxid") or "").strip(), f"manifest:{spp}")

            out_faa_raw = (r.get("out_faa") or "").strip()
            out_faa     = Path(out_faa_raw)

            if not out_faa.is_file():
                skipped_missing_file += 1
                warn(f"Skipping manifest row because out_faa is missing: {spp} -> {out_faa}")
                failed_rows.append({"species": spp, "reason": "out_faa_missing"})
                continue

            picked_assembly = (r.get(picked_col) or "").strip() or "NA"
            assembly_source = (r.get(source_col) or "").strip() or "NA"

            source_label = (r.get("source_label") or "").strip() or "REFSEQ"
            rows.append({
                "Species_name":    spp,
                "taxid":           taxid,
                "picked_assembly": picked_assembly,
                "assembly_source": assembly_source,
                "source_label":    source_label,
                "out_faa":         out_faa,
            })

    if not rows:
        # No usable OK rows is NOT a hard failure. A build can legitimately
        # have no Module I proteome source (e.g. a Module II / Module IV-only
        # database), which leaves this manifest header-only. It can also mean
        # every proteome row failed to resolve or download -- those failures
        # are already captured in failed_rows and surfaced in the build
        # summary. In both cases we warn and let the build proceed with
        # whatever other sources it has, instead of aborting the whole run.
        warn(
            f"No usable OK rows in {manifest_path} "
            f"(skipped_non_ok={skipped_non_ok}, skipped_missing_file={skipped_missing_file}); "
            f"continuing without proteome-derived proteins."
        )
        return rows, failed_rows

    info(
        f"Manifest loaded: {len(rows)} OK usable rows "
        f"(skipped_non_ok={skipped_non_ok}, skipped_missing_file={skipped_missing_file})"
    )

    return rows, failed_rows


# -----------------------------------------------------------
# USER-SUPPLIED METADATA LOADER
# -----------------------------------------------------------

USER_DEFAULTS = {
    "source":            "USER-supplied",
    "ncbi_taxid":        "32644",
    "genus":             "NA",
    "species":           "NA",
    "strain_or_gene":    "NA",
    "gene":              "NA",
    "protein_name":      "NA",
    "assembly":          "NA",
    "description":       "NA",
    "fetch_id":          "NA",
    "fetch_type":        "NA",
    "source_accession":  "NA",
    "returned_fasta_id": "NA",
}

REQUIRED_USER_COLS = {"protein_id"}


def load_user_metadata(path: Path):
    """
    Parse user_metadata.tsv.

    Required columns:
      protein_id

    Optional columns (defaults applied if missing or empty):
      source, ncbi_taxid, genus, species, strain_or_gene, gene, protein_name,
      assembly, description, fetch_id, fetch_type, source_accession,
      returned_fasta_id

    gene / protein_name are the explicit UniProt-style values exposed in the
    output header (GN= and the free-text protein name); both default to NA when
    not supplied, and may alternatively be carried on the FASTA header itself
    (GN= token / free-text description) for fetched proteins.
    """
    mapping = {}

    with path.open(newline="") as fh:
        # TSV metadata: tabs are the only delimiter and description fields hold
        # free text (brackets, semicolons, stray quotes from Excel). Disable
        # quote processing so a lone " never swallows subsequent rows.
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)

        if reader.fieldnames is None:
            raise RuntimeError(f"HARD FAIL: {path} has no header row.")

        missing_required = REQUIRED_USER_COLS - set(reader.fieldnames)

        if missing_required:
            raise RuntimeError(
                f"HARD FAIL: user metadata table is missing required columns: "
                f"{sorted(missing_required)}\nFound: {reader.fieldnames}"
            )

        for r in reader:
            pid = (r.get("protein_id") or "").strip()

            if not pid:
                warn("Skipping user metadata row with empty protein_id.")
                continue

            if pid in mapping:
                warn(f"Duplicate protein_id in user metadata, overwriting: {pid}")

            def get_field(col, _r=r):
                val = (_r.get(col) or "").strip()
                return val if val and val != "NA" else USER_DEFAULTS[col]

            raw_taxid = (r.get("ncbi_taxid") or "").strip()
            if raw_taxid and raw_taxid != "NA":
                taxid = require_taxid(raw_taxid, f"user_metadata:{pid}")
            else:
                taxid = USER_DEFAULTS["ncbi_taxid"]

            mapping[pid] = {
                "source":            get_field("source"),
                "ncbi_taxid":        taxid,
                "genus":             get_field("genus"),
                "species":           get_field("species"),
                "strain_or_gene":    get_field("strain_or_gene"),
                "gene":              get_field("gene"),
                "protein_name":      get_field("protein_name"),
                "assembly":          get_field("assembly"),
                "description":       get_field("description"),
                "fetch_id":          get_field("fetch_id"),
                "fetch_type":        get_field("fetch_type"),
                "source_accession":  get_field("source_accession"),
                "returned_fasta_id": get_field("returned_fasta_id"),
            }

    if not mapping:
        raise RuntimeError(
            f"HARD FAIL: no entries parsed from user metadata table: {path}"
        )

    info(f"User metadata loaded: {len(mapping)} entries from {path}")
    return mapping


# -----------------------------------------------------------
# NCBI EUTILS — SHARED HELPERS
# Mirrors the curl + XML parsing approach in Fetch_Species_INFO_2.0.sh
# Uses only Python standard library: urllib, xml.etree.ElementTree
# -----------------------------------------------------------

def _build_eutils_params(extra: dict, email: str = None, api_key: str = None) -> dict:
    params = {"tool": EUTILS_TOOL}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    params.update(extra)
    return params


def _eutils_get(endpoint: str, params: dict, retries: int = 5, sleep_sec: float = 2.0) -> str:
    """
    Fetch from NCBI eUtils and return the response body as a string.

    Retry behaviour mirrors fetch_to_file() in Fetch_Species_INFO_2.0.sh:
      - up to `retries` attempts
      - exponential back-off between attempts (2, 4, 8, 16 s)
      - returns None on total failure
    """
    url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                content = resp.read().decode("utf-8")
                if content:
                    return content
        except Exception as e:
            warn(f"eUtils fetch attempt {attempt}/{retries} failed ({endpoint}): {e}")

        if attempt < retries:
            time.sleep(2 ** attempt)

    return None


# -----------------------------------------------------------
# NCBI TAXONOMY LOOKUPS
# -----------------------------------------------------------

_NCBI_SLEEP_NO_KEY  = 2.0
_NCBI_SLEEP_API_KEY = 0.11


def _ncbi_sleep(api_key: str = None):
    time.sleep(_NCBI_SLEEP_API_KEY if api_key else _NCBI_SLEEP_NO_KEY)


def lookup_genus_species_from_taxid(
    taxid:     str,
    email:     str   = None,
    api_key:   str   = None,
    retries:   int   = 5,
    sleep_sec: float = 2.0,
) -> tuple:
    """
    Query NCBI taxonomy efetch to retrieve genus and species from a taxid.
    Mirrors get_taxonomy_scientific_name() in Fetch_Species_INFO_2.0.sh.
    Returns (genus, species) strings, or ("NA", "NA") on failure.
    """
    params  = _build_eutils_params(
        {"db": "taxonomy", "id": taxid, "retmode": "xml"},
        email=email, api_key=api_key,
    )
    content = _eutils_get("efetch.fcgi", params, retries=retries, sleep_sec=sleep_sec)

    if not content:
        warn(f"NCBI efetch returned no content for taxid {taxid}.")
        return "NA", "NA"

    try:
        root = ET.fromstring(content)
    except ET.ParseError as e:
        warn(f"XML parse error for taxid {taxid}: {e}")
        return "NA", "NA"

    taxon = root.find("Taxon")
    if taxon is None:
        warn(f"No <Taxon> element in NCBI response for taxid {taxid}.")
        return "NA", "NA"

    rank     = (taxon.findtext("Rank")           or "").strip()
    sci_name = (taxon.findtext("ScientificName") or "").strip()

    genus   = "NA"
    species = "NA"

    lineage_ex = taxon.find("LineageEx")
    if lineage_ex is not None:
        for node in lineage_ex.findall("Taxon"):
            node_rank = (node.findtext("Rank")           or "").strip()
            node_name = (node.findtext("ScientificName") or "").strip()

            if node_rank == "genus":
                genus = node_name
            elif node_rank == "species":
                parts   = node_name.split()
                species = " ".join(parts[1:]) if len(parts) > 1 else node_name

    if rank == "species" and sci_name:
        parts = sci_name.split()
        if genus == "NA" and parts:
            genus = parts[0]
        if species == "NA":
            species = " ".join(parts[1:]) if len(parts) > 1 else sci_name
    elif rank == "genus" and sci_name:
        if genus == "NA":
            genus = sci_name

    if genus == "NA":
        warn(
            f"Could not resolve genus from NCBI for taxid {taxid} "
            f"(rank={rank}, ScientificName={sci_name})."
        )

    return genus, species


def lookup_taxid_from_name(
    organism_name: str,
    email:     str   = None,
    api_key:   str   = None,
    retries:   int   = 5,
    sleep_sec: float = 2.0,
) -> str:
    """
    Query NCBI taxonomy esearch to find taxid from an organism name.
    Uses the same query format as get_taxid_from_taxonomy() in
    Fetch_Species_INFO_2.0.sh:
      "{name}"[Scientific Name] OR "{name}"[All Names]
    Returns taxid string, or None on failure/no match.
    """
    query = (
        f'"{organism_name}"[Scientific Name] '
        f'OR "{organism_name}"[All Names]'
    )

    params  = _build_eutils_params(
        {"db": "taxonomy", "term": query, "retmode": "xml", "retmax": "1"},
        email=email, api_key=api_key,
    )
    content = _eutils_get("esearch.fcgi", params, retries=retries, sleep_sec=sleep_sec)

    if not content:
        warn(f"NCBI esearch returned no content for organism '{organism_name}'.")
        return None

    try:
        root = ET.fromstring(content)
    except ET.ParseError as e:
        warn(f"XML parse error for organism search '{organism_name}': {e}")
        return None

    id_list = root.find("IdList")
    if id_list is None:
        warn(f"No <IdList> in NCBI esearch response for '{organism_name}'.")
        return None

    ids = [el.text.strip() for el in id_list.findall("Id") if el.text]

    if not ids:
        warn(f"No NCBI taxid found for organism name: '{organism_name}'.")
        return None

    if len(ids) > 1:
        warn(
            f"Multiple taxids found for '{organism_name}': {ids}. "
            f"Using first hit ({ids[0]}). "
            f"Provide a more specific name or supply taxid directly for precision."
        )

    return ids[0]


def enrich_user_metadata(
    user_meta:  dict,
    ncbi_email: str,
    api_key:    str   = None,
    sleep_sec:  float = 2.0,
) -> dict:
    """
    For each metadata row, cross-fill missing taxonomy fields using NCBI eUtils.

    Rules (mirrors PATH A / PATH B logic in Fetch_Species_INFO_2.0.sh):
      taxid known,   genus/species missing -> efetch  genus + species from taxid
      taxid unknown, genus present         -> esearch taxid from genus [+ species]
      both present                         -> use as-is, no lookup
      neither present                      -> warn only, keep defaults

    Uses in-memory caches so repeated taxids / organism names are only queried
    once per build. Failed lookups are cached too, which prevents repeated
    failed NCBI calls for the same unresolved organism.

    Modifies user_meta in-place and returns it.
    """
    warn(
        f"Starting NCBI taxonomy enrichment for {len(user_meta)} user metadata rows "
        f"(sleep={sleep_sec}s between requests, "
        f"{'API key set' if api_key else 'no API key'})."
    )

    # Caches are per run/build. They do not alter user-provided metadata files.
    taxid_to_org_cache = {}   # taxid string -> (genus, species)
    name_to_taxid_cache = {}  # normalized organism query -> taxid string or ""

    enriched         = 0
    lookup_failed    = 0
    skipped_complete = 0
    skipped_empty    = 0
    ncbi_queries     = 0
    cache_hits_taxid = 0
    cache_hits_name  = 0

    for pid, m in user_meta.items():
        taxid   = m["ncbi_taxid"]
        genus   = m["genus"]
        species = m["species"]

        taxid_known   = taxid != USER_DEFAULTS["ncbi_taxid"]
        genus_known   = genus != "NA"
        species_known = species != "NA"
        org_complete  = genus_known and species_known

        if taxid_known and not org_complete:
            used_cache = False

            if taxid in taxid_to_org_cache:
                g, s = taxid_to_org_cache[taxid]
                cache_hits_taxid += 1
                used_cache = True
            else:
                warn(
                    f"[{pid}] taxid={taxid} present but genus/species incomplete; "
                    f"querying NCBI efetch..."
                )
                g, s = lookup_genus_species_from_taxid(
                    taxid, email=ncbi_email, api_key=api_key, sleep_sec=sleep_sec,
                )
                taxid_to_org_cache[taxid] = (g, s)
                ncbi_queries += 1
                _ncbi_sleep(api_key)

            if g != "NA":
                if not genus_known:
                    m["genus"] = g
                if not species_known:
                    m["species"] = s
                enriched += 1
                if not used_cache:
                    warn(f"  -> genus={m['genus']}, species={m['species']}")
            else:
                lookup_failed += 1
                if not used_cache:
                    warn(f"  -> NCBI lookup failed; genus/species remain NA for '{pid}'.")

        elif not taxid_known and genus_known:
            query = f"{genus} {species}".strip() if species_known else genus
            cache_key = " ".join(query.lower().split())
            used_cache = False

            if cache_key in name_to_taxid_cache:
                found_taxid = name_to_taxid_cache[cache_key]
                cache_hits_name += 1
                used_cache = True
            else:
                warn(
                    f"[{pid}] organism='{query}' present but taxid unknown; "
                    f"querying NCBI esearch..."
                )
                found_taxid = lookup_taxid_from_name(
                    query, email=ncbi_email, api_key=api_key, sleep_sec=sleep_sec,
                )
                # Store "" for failed lookups so repeated failures are not re-queried.
                name_to_taxid_cache[cache_key] = found_taxid or ""
                ncbi_queries += 1
                _ncbi_sleep(api_key)

            if found_taxid:
                m["ncbi_taxid"] = found_taxid
                enriched += 1
                if not used_cache:
                    warn(f"  -> taxid={found_taxid}")
            else:
                lookup_failed += 1
                if not used_cache:
                    warn(
                        f"  -> NCBI lookup failed; taxid remains "
                        f"{USER_DEFAULTS['ncbi_taxid']} for '{pid}'."
                    )

        elif taxid_known and org_complete:
            skipped_complete += 1

        else:
            skipped_empty += 1
            warn(
                f"[{pid}] No taxid or organism supplied; all taxonomy fields "
                f"remain at defaults. Consider adding values to the metadata table."
            )

    warn(
        f"NCBI enrichment complete: "
        f"enriched={enriched}, "
        f"lookup_failed={lookup_failed}, "
        f"already_complete={skipped_complete}, "
        f"no_info_supplied={skipped_empty}, "
        f"ncbi_queries={ncbi_queries}, "
        f"cache_hits_taxid={cache_hits_taxid}, "
        f"cache_hits_name={cache_hits_name}."
    )

    return user_meta


# -----------------------------------------------------------
# FASTA HEADER PROVENANCE HELPERS
# -----------------------------------------------------------

def parse_fasta_header_keyvals(header: str) -> dict:
    """Parse whitespace-delimited key=value tokens from a FASTA header."""
    kv = {}
    for token in (header or "").split()[1:]:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key:
            kv[key] = value
    return kv


_ORG_BRACKET_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


def harvest_organism_from_desc(desc: str):
    """Pull the trailing [Organism] from an NCBI-style protein title."""
    if not desc:
        return "", "", ""
    m = _ORG_BRACKET_RE.search(desc.strip())
    if not m:
        return "", "", ""
    organism = m.group(1).strip()
    toks = organism.split()
    genus = toks[0] if toks else ""
    species = toks[1] if len(toks) > 1 else ""
    return genus, species, organism


def clean_product_desc(fasta_desc: str) -> str:
    """Drop leading key=value provenance tokens, keep the human product text."""
    toks = (fasta_desc or "").split()
    return " ".join(t for t in toks if "=" not in t).strip()


def clean_protein_name(text: str) -> str:
    """Human-readable product/protein name, with the UniProt key=value tail
    (OS=/OX=/GN=/PE=/SV=...), any trailing [Organism], and stray key=value
    provenance tokens removed. Returns 'NA' when nothing meaningful is left.

    Works for every source's raw description:
      UniProt : 'Catalase OS=Homo sapiens OX=9606 GN=CAT ...' -> 'Catalase'
      NCBI faa: 'catalase [Homo sapiens]'                     -> 'catalase'
      PROKKA  : 'hypothetical protein'                        -> 'hypothetical protein'
    """
    if not text:
        return "NA"
    name = re.split(r"(?:^|\s)[A-Z]{2}=", text)[0]   # cut UniProt key tail
    # Drop ALL trailing [..] groups (NCBI [Organism], HOMD's extra [HMT-579], etc.)
    prev = None
    while prev != name:
        prev = name
        name = _ORG_BRACKET_RE.sub("", name).rstrip()
    name = " ".join(t for t in name.split() if "=" not in t)
    name = safe_one_line(name)
    return name or "NA"


def seed_fetched_protein_stubs(user_meta: dict, fasta_path) -> None:
    """Auto-build metadata for fetched proteins straight from the fetched FASTA.

    For each fetched record (defline carries fetch provenance) this ensures a
    user_meta row exists and backfills genus/species/description/source from
    the returned NCBI title and header provenance where they are still NA.
    Fallback-only: never overwrites a value the user supplied.
    """
    seeded = harvested = 0
    for header, _seq in iter_fasta(fasta_path):
        parts = header.split()
        if not parts:
            continue
        pid = parts[0]
        kv = parse_fasta_header_keyvals(header)

        if not (kv.get("returned_fasta_id") or kv.get("fetch_id") or kv.get("fetch_type")):
            continue

        desc = " ".join(parts[1:]) if len(parts) > 1 else ""
        product = clean_product_desc(desc)
        genus, species, _organism = harvest_organism_from_desc(product)

        m = user_meta.get(pid)
        if m is None:
            m = dict(USER_DEFAULTS)
            m["fetch_id"]          = kv.get("fetch_id", "NA") or "NA"
            m["fetch_type"]        = kv.get("fetch_type", "NA") or "NA"
            m["source_accession"]  = kv.get("source_accession", pid) or pid
            m["returned_fasta_id"] = kv.get("returned_fasta_id", pid) or pid
            # Prefer the registry source_label carried on the fetched header
            # (written by fetchProteinAccessions_batch.py) over the generic
            # USER_DEFAULTS["source"] fallback, so metadata-free protein
            # sources still get correctly attributed in the manifest instead
            # of collapsing into a shared "USER-supplied" bucket.
            if kv.get("source_label"):
                m["source"] = kv["source_label"]
            user_meta[pid] = m
            seeded += 1

        if m.get("genus", "NA") == "NA" and genus:
            m["genus"] = genus
            harvested += 1
        if m.get("species", "NA") == "NA" and species:
            m["species"] = species
        if m.get("description", "NA") == "NA" and product:
            m["description"] = product

    if seeded or harvested:
        info(f"Fetched-protein autofill: stubbed {seeded} rows, "
             f"harvested organism for {harvested}.")


def choose_provenance_value(*values, default="NA"):
    """Return the first non-empty/non-NA provenance value."""
    for value in values:
        value = (value or "").strip()
        if value and value != "NA":
            return value
    return default


# -----------------------------------------------------------
# WRITERS
# -----------------------------------------------------------

# Canonical manifest column order.  All 14 fields default to "NA" so each
# call site only needs to supply the values it actually has.
MANIFEST_COLUMNS = [
    "seqid",
    "source",
    "protein_id",
    "fetch_id",
    "fetch_type",
    "source_accession",
    "returned_fasta_id",
    "assembly_or_uniprot",
    "hmt_id",
    "ncbi_taxid",
    "genus",
    "species",
    "strain_or_gene",
    "gene_name",
    "protein_name",
    "description",
    "lineage_taxid",
]

# Sources whose proteins should be tallied under a different taxonomic
# lineage than their literal NCBI-resolved organism, for lineage-based
# summary/rollup purposes only. For example, HERV loci correctly resolve
# to their literal GenBank host organism "Homo sapiens" (they genuinely
# are annotated under that organism in GenBank) but should be counted as
# "Human endogenous retroviruses" rather than human in any species/genus
# rollup. ncbi_taxid/genus/species are NEVER modified by this -- they stay
# the true NCBI-resolved values for provenance/audit. This dict only
# affects the separate lineage_taxid column.
#
# Populated entirely from the registry -- ANY source can declare
# "lineage_taxid_override=<taxid>" in its options column (same syntax as
# fetch_source=) and it lands here automatically via load_lineage_overrides()
# below. No code change is needed to add a new lineage-override source.
#
# Applied centrally in write_manifest_row() (keyed on row.source, the
# literal source_label) rather than at each ManifestRow construction site,
# so no current or future fetch path can forget to apply it.
SOURCE_LINEAGE_TAXID_OVERRIDES: dict[str, str] = {}


def load_lineage_overrides(source_plan_path) -> dict[str, str]:
    """Read lineage_taxid_override=<taxid> options from the normalized
    source plan, keyed by source_label. Lets any registry source declare a
    file-wide lineage override via its own options column, e.g.:

        options = fetch_source=ncbi;lineage_taxid_override=206037

    Silently returns {} if no source plan was supplied -- the build still
    runs fine, it just has no lineage overrides (every row's lineage_taxid
    falls back to its true ncbi_taxid, same as if this feature didn't exist).
    """
    overrides: dict[str, str] = {}
    if not source_plan_path:
        return overrides
    source_plan_path = Path(source_plan_path)
    if not source_plan_path.is_file():
        warn(f"--source_plan given but not found: {source_plan_path}; "
             f"no lineage overrides will be applied.")
        return overrides
    with source_plan_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            label = (row.get("source_label") or "").strip()
            opts_raw = (row.get("options") or "").strip()
            if not label or not opts_raw:
                continue
            for item in opts_raw.split(";"):
                item = item.strip()
                if not item or "=" not in item:
                    continue
                k, v = item.split("=", 1)
                if k.strip().lower() == "lineage_taxid_override" and v.strip():
                    overrides[label] = v.strip()
                    info(f"Lineage override from registry: source={label!r} "
                         f"-> lineage_taxid={v.strip()!r}")
    return overrides


@dataclass
class ManifestRow:
    seqid:              str = "NA"
    source:             str = "NA"
    protein_id:         str = "NA"
    fetch_id:           str = "NA"
    fetch_type:         str = "NA"
    source_accession:   str = "NA"
    returned_fasta_id:  str = "NA"
    assembly_or_uniprot: str = "NA"
    hmt_id:             str = "NA"
    ncbi_taxid:         str = "NA"
    genus:              str = "NA"
    species:            str = "NA"
    strain_or_gene:     str = "NA"
    gene_name:          str = "NA"
    protein_name:       str = "NA"
    description:        str = "NA"
    lineage_taxid:       str = "NA"


def next_internal_id(prefix: str, uid: int):
    return f"{prefix}-{uid:011d}"


def write_manifest_header(mout) -> None:
    """Write the manifest TSV header row."""
    mout.write("\t".join(MANIFEST_COLUMNS) + "\n")


def write_manifest_row(mout, row: ManifestRow) -> None:
    """Write one manifest row in canonical column order.

    lineage_taxid defaults to the row's true ncbi_taxid; if the row's
    source has a registry-declared lineage override
    (SOURCE_LINEAGE_TAXID_OVERRIDES, populated from lineage_taxid_override=
    in that source's options column) that is applied instead.
    ncbi_taxid/genus/species are left untouched -- only this separate
    column is affected.
    """
    if row.lineage_taxid in ("NA", "", None):
        row.lineage_taxid = SOURCE_LINEAGE_TAXID_OVERRIDES.get(row.source, row.ncbi_taxid)
    mout.write("\t".join(str(getattr(row, c)) for c in MANIFEST_COLUMNS) + "\n")


def write_fasta_record(fout, header_line: str, seq: str):
    fout.write(header_line + "\n")
    for chunk in wrap_sequence(seq):
        fout.write(chunk + "\n")


# -----------------------------------------------------------
# UNIFIED, PIPE-FREE FASTA DEFLINE
# -----------------------------------------------------------
#
# Every source funnels through format_defline() so the header format can never
# diverge per-source again. The layout is deliberately PIPE-FREE: the unique
# seqid is the sole first-token identifier and all metadata is space-delimited
# key=value, with UniProt-recognised tags (OS=/OX=/GN=) first so MaxQuant and
# FragPipe/Philosopher read organism + gene, followed by maniFasta provenance
# tags. No '|' is ever emitted, so the default UniProt pipe-splitting parse
# rules can't mis-assign protein IDs.
#
#   >{seqid} {protein_name} OS={organism} OX={taxid} GN={gene} \
#    src_db={source} group={group} version={version} acc={accession} \
#    assembly={assembly} src={asm_src}[ HMT={hmt}]

# Central registry provenance, populated once from the normalized source plan
# (same pattern as SOURCE_LINEAGE_TAXID_OVERRIDES). Keyed by source_label so any
# lane can stamp the right group/version without threading values per-fetcher.
SOURCE_REGISTRY_META: dict[str, dict] = {}     # source_label -> {"group","version"}
COLLECTION_LABELS:    dict[str, str]  = {}      # COLLECTION(upper) -> source_label


def load_source_registry_meta(source_plan_path) -> dict[str, dict]:
    """Read source_label -> {group, version} (and, as a side effect, the
    collection -> source_label map for supported_collection lanes) from the
    normalized source plan. Silently returns {} if no plan is supplied -- the
    build still runs, headers just carry group=NA / version=NA."""
    meta: dict[str, dict] = {}
    if not source_plan_path:
        return meta
    source_plan_path = Path(source_plan_path)
    if not source_plan_path.is_file():
        warn(f"--source_plan given but not found: {source_plan_path}; "
             f"source_group/version will be NA in headers.")
        return meta
    with source_plan_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            label = (row.get("source_label") or "").strip()
            if label:
                meta[label] = {
                    "group":   (row.get("source_group") or "").strip() or "NA",
                    "version": (row.get("version") or "").strip() or "NA",
                }
                coll = (row.get("collection") or "").strip()
                if coll and (row.get("enabled") or "").strip().upper() == "TRUE":
                    # case-insensitive key so registry 'cRAP' resolves under 'CRAP'
                    COLLECTION_LABELS[coll.upper()] = label
    info(f"Source registry meta loaded for {len(meta)} labels "
         f"({len(COLLECTION_LABELS)} collection labels).")
    return meta


def _hdr_safe(v) -> str:
    """Header-token value: collapse internal whitespace and strip the pipe so a
    '|' can NEVER reach an output defline (downstream tools split UniProt-style
    headers on '|'). Pipes only ever appear in external inputs (UniProt
    sp|acc|name, RCSB, NCBI gi|...); this is the last line of defence."""
    v = safe_one_line(str(v)) if v is not None else ""
    return v.replace("|", "_")


def format_defline(
    seqid: str,
    source: str,
    accession: str,
    organism: str,
    taxid: str,
    *,
    assembly: str = "NA",
    asm_src: str = "NA",
    gene: str = "NA",
    protein_name: str = "NA",
    hmt: str = "NA",
    source_group: str = None,   # None -> resolve from registry by source_label
    version: str = None,        # None -> resolve from registry by source_label
) -> str:
    """THE maniFasta defline builder (pipe-free, UniProt-style key=value)."""

    def kv(v) -> str:
        v = _hdr_safe(v)
        return v if (v and v != "NA") else "NA"

    # Resolve registry provenance centrally unless the caller overrode it.
    meta = SOURCE_REGISTRY_META.get(source, {})
    if source_group is None:
        source_group = meta.get("group", "NA")
    if version is None:
        version = meta.get("version", "NA")

    tags = [
        f"OS={_hdr_safe(organism)}",
        f"OX={kv(taxid)}",
        f"GN={kv(gene)}",
        f"src_db={kv(source)}",
        f"group={kv(source_group)}",
        f"version={kv(version)}",
        f"acc={kv(accession)}",
        f"assembly={kv(assembly)}",
        f"src={kv(asm_src)}",
    ]
    if kv(hmt) != "NA":
        tags.append(f"HMT={kv(hmt)}")

    parts = [f">{_hdr_safe(seqid)}"]
    name = _hdr_safe(protein_name)
    if name and name != "NA":
        parts.append(name)        # free text in UniProt's slot, omitted if NA
    parts.extend(tags)
    return " ".join(parts)


# -----------------------------------------------------------
# MAIN DATABASE BUILD
# -----------------------------------------------------------


# -----------------------------------------------------------
# BUILD SUMMARY
# -----------------------------------------------------------

def write_build_summary(path: Path, args, source_counts: dict,
                        total: int, start_time, run_label: str,
                        failed_sources: list | None = None):
    """Write a human-readable TSV summary of the build."""
    import datetime
    end_time  = datetime.datetime.now()
    duration  = str(end_time - start_time).split(".")[0]  # HH:MM:SS
    timestamp = start_time.strftime("%Y-%m-%d %H:%M:%S")

    with path.open("w", encoding="utf-8") as f:
        def row(key, val):
            f.write(f"{key}\t{val}\n")

        f.write("# maniFasta build summary\n")
        f.write(f"# Generated: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("#\n")
        f.write("field\tvalue\n")

        # Run metadata
        row("run_label",        run_label or "unset")
        row("run_timestamp",    timestamp)
        row("run_duration",     duration)
        row("total_proteins",   total)

        # Surface the missing/skipped count up top alongside the protein total
        # so it is visible without scrolling to the failed-sources section.
        fs_list = failed_sources or []
        row("total_missing_accessions", len(fs_list))

        f.write("#\n# ── Sources ──\n")
        # Fixed sources
        if args.homd_fasta:
            row("source_HOMD",   source_counts.get("HOMD", 0))
        if args.crap_fasta:
            crap_label = Path(args.crap_fasta).stem
            row(f"source_CRAP_{crap_label}", source_counts.get("CRAP", 0))

        # User-defined sources from manifest (by source_label)
        manifest_sources = {k: v for k, v in source_counts.items()
                            if k not in ("HOMD", "CRAP",
                                         "_manifest_total", "_user_total")}
        for label, count in sorted(manifest_sources.items()):
            row(f"source_{label}", count)

        # User-supplied subtotal (mod_III sources, already itemised above)
        if "_user_total" in source_counts:
            row("subtotal_user_supplied", source_counts["_user_total"])

        f.write("#\n# ── Inputs ──\n")
        row("homd_fasta",          str(args.homd_fasta)  if args.homd_fasta  else "not included")
        row("homd_tax",            str(args.homd_tax)    if args.homd_tax    else "not included")
        row("gca_info",            str(args.gca_info)    if args.gca_info    else "not included")
        row("human_fasta",         str(args.human_fasta) if args.human_fasta else "not included")
        row("crap_fasta",          str(args.crap_fasta)  if args.crap_fasta  else "not included")
        row("proteomes_manifest",  str(args.proteomes_manifest) if args.proteomes_manifest else "not included")
        row("user_fasta",          str(args.user_fasta)  if args.user_fasta  else "not included")
        row("user_metadata",       str(args.user_metadata) if args.user_metadata else "not included")

        f.write("#\n# ── Failed / skipped sources ──\n")
        row("failed_source_count", len(fs_list))
        for fs in fs_list:
            row(f"FAILED\t{fs['species']}", fs["reason"])

        f.write("#\n# ── Outputs ──\n")
        row("out_fasta",           str(args.out_fasta))
        row("out_map",             str(args.out_map))

    info(f"Build summary written: {path}")


def build_db(args):
    uid = 1
    import datetime
    start_time = datetime.datetime.now()

    SOURCE_LINEAGE_TAXID_OVERRIDES.update(
        load_lineage_overrides(getattr(args, "source_plan", None))
    )
    # Populates SOURCE_REGISTRY_META (group/version) and COLLECTION_LABELS.
    SOURCE_REGISTRY_META.update(
        load_source_registry_meta(getattr(args, "source_plan", None))
    )
    source_counts = {}   # source_label -> protein count
    failed_sources: list = []  # [{species, reason}] for build summary

    # Fold in accessions a fetcher skipped under --allow-missing so they appear
    # in the build summary's failed/skipped section. skipped_ids holds the bare
    # accessions so they can be excused from the user-metadata "unused rows"
    # hard-fail below (their metadata row exists but the FASTA record does not).
    skipped_ids: set = set()
    if getattr(args, "failed_accessions", None) and Path(args.failed_accessions).is_file():
        with Path(args.failed_accessions).open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                label  = (r.get("source_label") or "NA").strip()
                acc    = (r.get("accession") or "NA").strip()
                reason = (r.get("reason") or "missing").strip()
                failed_sources.append({"species": f"{label}:{acc}", "reason": reason})
                if acc and acc != "NA":
                    skipped_ids.add(acc)
        if skipped_ids:
            warn(f"Loaded {len(skipped_ids)} skipped/missing accession(s) for the build summary.")

    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.out_map.parent.mkdir(parents=True, exist_ok=True)

    with args.out_fasta.open("w") as fout, args.out_map.open("w") as mout:

        write_manifest_header(mout)

        # ---------------------------------------------------
        # HOMD (optional)
        #
        # All three HOMD args must be supplied together or not at all.
        # If a body-site-filtered taxonomy TSV is passed as --homd_tax,
        # proteins whose HMT ID is absent from the filtered table are
        # silently skipped — this IS the filtering mechanism.
        # ---------------------------------------------------
        homd_args_provided = (
            args.homd_fasta is not None
            and args.gca_info  is not None
            and args.homd_tax  is not None
        )
        homd_args_partial = (
            (args.homd_fasta is not None)
            + (args.gca_info  is not None)
            + (args.homd_tax  is not None)
        )

        if 0 < homd_args_partial < 3:
            raise RuntimeError(
                "HARD FAIL: --homd_fasta, --gca_info, and --homd_tax must all be "
                "supplied together, or all omitted. "
                f"Got {homd_args_partial}/3 arguments."
            )

        if homd_args_provided:
            gca_to_hmt = load_gca_to_hmt(args.gca_info)
            hmt_to_tax = load_hmt_to_tax(args.homd_tax)

            info(
                f"HOMD taxonomy table loaded: {len(hmt_to_tax)} HMT entries "
                f"(from {args.homd_tax})."
            )

            homd_written        = 0
            homd_filtered_out   = 0   # skipped because HMT not in (filtered) tax table
            homd_no_assembly    = 0   # HARD FAIL counter — kept for diagnostic clarity

            for header, seq in iter_fasta(args.homd_fasta):
                parts = header.split()

                if not parts:
                    continue

                prot_id  = parts[0]
                # PROKKA locus tags contain underscores, so anchor on the leading
                # GC[AF]_<digits>.<ver> rather than rsplit('_', 1).
                m_asm    = re.match(r"^(GC[AF]_\d+(?:\.\d+)?)_", prot_id)
                assembly = m_asm.group(1) if m_asm else prot_id.rsplit("_", 1)[0]

                # Missing assembly is always a data integrity error
                if assembly not in gca_to_hmt:
                    raise RuntimeError(
                        f"HARD FAIL: HOMD assembly '{assembly}' (from protein "
                        f"'{prot_id}') not found in GCA_ID_info file: {args.gca_info}"
                    )

                gmeta       = gca_to_hmt[assembly]
                hmt         = norm_hmt(gmeta["hmt"])
                genus_gca   = gmeta["genus"]
                species_gca = gmeta["species"]
                strain      = gmeta["strain"]

                # HMT absent from taxonomy table = filtered out by homd_filter.py
                # (or a genuine data gap). Skip with a warning rather than hard-failing
                # so body-site filtering works correctly.
                hmt_key = norm_hmt(hmt)

                if hmt_key not in hmt_to_tax:
                    homd_filtered_out += 1
                    continue

                tmeta    = hmt_to_tax[hmt_key]
                taxid    = require_taxid(tmeta["taxid"], f"HOMD:{prot_id}")

                genus    = (tmeta["genus"]   or genus_gca).strip()   or genus_gca
                species  = (tmeta["species"] or species_gca).strip() or species_gca
                organism = require_organism(f"{genus} {species}", f"HOMD:{prot_id}")

                internal_id = next_internal_id(DB_PREFIX, uid)
                uid += 1

                desc         = " ".join(parts[1:]) if len(parts) > 1 else ""
                # Gene symbol injected as GN=<gene> by the HOMD fetcher from the
                # per-genome PROKKA .tsv (locus_tag -> gene). The .faa itself
                # carries only the product, so absent GN= -> gene stays NA.
                m_gene       = re.search(r"(?:^|\s)GN=(\S+)", desc)
                gene         = m_gene.group(1) if m_gene else "NA"
                protein_name = clean_protein_name(desc)   # also strips the GN= tail

                # GCA CONSOLIDATION: HOMD PROKKA locus tags already carry the
                # assembly (GCA_x_00001); the fetcher historically prepended the
                # GCA again, yielding GCA_x_GCA_x_00001. Collapse any repeated
                # assembly prefix down to ONE so protein_id == GCA_x_00001, and
                # carry the bare assembly only in the single assembly_or_uniprot
                # column. Works whether the input is singly or doubly prefixed.
                locus = prot_id
                dup = assembly + "_" + assembly + "_"
                if locus.startswith(dup):
                    locus = assembly + "_" + locus[len(dup):]   # GCA_x_GCA_x_00001 -> GCA_x_00001

                # Registry label for the HOMD collection (falls back to "HOMD").
                homd_label = COLLECTION_LABELS.get("HOMD", "HOMD")

                out_header = format_defline(
                    internal_id, homd_label, locus, organism, taxid,
                    assembly=assembly, asm_src="homd_prokka",
                    gene=gene, protein_name=protein_name, hmt=hmt,
                )

                write_fasta_record(fout, out_header, seq)

                write_manifest_row(mout, ManifestRow(
                    seqid=internal_id,
                    source=homd_label,
                    protein_id=locus,                 # GCA_x_00001 (single GCA)
                    fetch_type="HOMD",
                    source_accession=assembly,        # the consolidated GCA column
                    returned_fasta_id=locus,
                    assembly_or_uniprot=assembly,
                    hmt_id=hmt,
                    ncbi_taxid=taxid,
                    genus=genus,
                    species=species,
                    strain_or_gene=strain,
                    gene_name=gene,                   # from PROKKA .tsv (GN=), else NA
                    protein_name=protein_name,
                    description=protein_name,
                ))

                homd_written += 1

            if homd_filtered_out > 0:
                warn(
                    f"HOMD proteins skipped (HMT not in taxonomy table / filtered out): "
                    f"{homd_filtered_out}"
                )

            if homd_written == 0:
                raise RuntimeError(
                    "HARD FAIL: HOMD args were supplied but 0 proteins were written. "
                    "If a body-site filter was applied, check that the requested "
                    "site(s) produced at least one matching HMT entry."
                )

            info(
                f"HOMD proteins written: {homd_written} "
                f"(filtered out: {homd_filtered_out})"
            )
            source_counts["HOMD"] = homd_written

        else:
            info("HOMD not included in this build (no HOMD arguments supplied).")

        # ---------------------------------------------------
        # HUMAN / SARS2 / CRAP
        # ---------------------------------------------------

        # Known taxonomy for cRAP records that lack OS=/OX= in their headers.
        # Keys are the protein id token (second pipe-delimited field).
        # Sources:
        #   cRAP126: Endoproteinase GluC = Staphylococcal V8 Protease (S. aureus)
        #            https://www.neb.com/en-us/products/p8100-endoproteinase-gluc
        #   cRAP127: rLys-C = Protease IV from Pseudomonas aeruginosa
        #            https://www.promega.com/products/mass-spectrometry/proteases-and-surfactants/rlys-c-mass-spec-grade/
        CRAP_TAXON_OVERRIDES: dict[str, tuple[str, str]] = {
            "cRAP126": ("1280", "Staphylococcus aureus"),
            "cRAP127": ("287",  "Pseudomonas aeruginosa"),
        }

        def process_uniprot_like(source_label: str, fasta_path: Path):
            nonlocal uid

            written = 0
            skipped = 0

            for header, seq in iter_fasta(fasta_path):
                try:
                    p = parse_uniprot_like_header(header)
                except RuntimeError as e:
                    if source_label.upper() == "CRAP":
                        # Extract the protein id token (sp|cRAP126|... -> cRAP126)
                        parts = header.lstrip(">").split("|")
                        prot_token = parts[1].strip() if len(parts) >= 2 else parts[0].split()[0]
                        if prot_token in CRAP_TAXON_OVERRIDES:
                            override_taxid, override_organism = CRAP_TAXON_OVERRIDES[prot_token]
                            desc = " ".join(header.split()[1:]) if len(header.split()) > 1 else ""
                            p = {
                                "id":          prot_token,
                                "gene":        "NA",
                                "taxid":       override_taxid,
                                "organism":    override_organism,
                                "description": desc,
                            }
                            info(
                                f"cRAP record '{prot_token}': applied taxon override "
                                f"({override_organism}, taxid={override_taxid})"
                            )
                        else:
                            skipped += 1
                            warn(
                                f"Skipping CRAP record due to missing OS/OX: "
                                f"{e} | header='{header[:120]}'"
                            )
                            continue
                    else:
                        raise

                prot     = p["id"]
                gene     = p["gene"]
                taxid    = require_taxid(p["taxid"],       f"{source_label}:{prot}")
                organism = require_organism(p["organism"], f"{source_label}:{prot}")
                desc     = p["description"]

                genus, species_only = split_organism_for_map(
                    organism, fallback_genus=source_label,
                )
                protein_name = clean_protein_name(desc)
                internal_id = next_internal_id(DB_PREFIX, uid)
                uid += 1

                out_header = format_defline(
                    internal_id, source_label, prot, organism, taxid,
                    assembly="NA", asm_src="NA",      # accession-level, no assembly
                    gene=gene, protein_name=protein_name,
                )

                write_fasta_record(fout, out_header, seq)

                write_manifest_row(mout, ManifestRow(
                    seqid=internal_id,
                    source=source_label,
                    protein_id=prot,
                    fetch_type=source_label,
                    source_accession=prot,
                    returned_fasta_id=prot,
                    assembly_or_uniprot=prot,
                    ncbi_taxid=taxid,
                    genus=genus,
                    species=species_only,
                    strain_or_gene=gene,
                    gene_name=gene,
                    protein_name=protein_name,
                    description=protein_name,
                ))

                written += 1

            if written == 0:
                raise RuntimeError(
                    f"HARD FAIL: wrote 0 proteins from {source_label} FASTA: {fasta_path}"
                )

            info(f"{source_label} proteins written: {written}; skipped: {skipped}")
            source_counts[source_label] = source_counts.get(source_label, 0) + written

        if args.human_fasta is not None:
            process_uniprot_like(COLLECTION_LABELS.get("HUMAN_UNIPROT", "HUMAN"), args.human_fasta)

        if args.crap_fasta is not None:
            process_uniprot_like(COLLECTION_LABELS.get("CRAP", "cRAP"), args.crap_fasta)

        # ---------------------------------------------------
        # EXTRA PROTEOMES FROM MANIFEST
        # ---------------------------------------------------
        if args.proteomes_manifest is not None:
            manifest_rows, manifest_failed = load_proteomes_manifest(args.proteomes_manifest)
            failed_sources.extend(manifest_failed)
            refseq_written = 0

            for r in manifest_rows:
                spp       = r["Species_name"]
                taxid     = r["taxid"]
                asm       = r["picked_assembly"] or "NA"
                asm_src   = r["assembly_source"] or "NA"
                src_label = r.get("source_label") or "REFSEQ"
                faa       = r["out_faa"]

                organism = require_organism(spp, f"manifest_species:{asm}:{faa.name}")

                for header, seq in iter_fasta(faa):
                    parts   = header.split()
                    prot_id = parts[0] if parts else "NA"

                    if prot_id == "NA":
                        raise RuntimeError(
                            f"HARD FAIL: could not parse protein id in {faa} "
                            f"(header='{header[:80]}...')"
                        )

                    desc         = " ".join(parts[1:]) if len(parts) > 1 else ""
                    protein_name = clean_protein_name(desc)
                    internal_id  = next_internal_id(DB_PREFIX, uid)
                    uid += 1

                    # NCBI deflines rarely carry a gene symbol; harvest [gene=XXX]
                    # when present, else GN stays NA (honest "where possible").
                    m_gn = re.search(r"\[gene=([^\]]+)\]", desc)
                    gene = m_gn.group(1).strip() if m_gn else "NA"

                    out_header = format_defline(
                        internal_id, src_label, prot_id, organism, taxid,
                        assembly=asm, asm_src=asm_src,
                        gene=gene, protein_name=protein_name,
                    )

                    write_fasta_record(fout, out_header, seq)

                    genus, species_only = split_organism_for_map(
                        organism, fallback_genus=src_label,
                    )

                    source_counts[src_label] = source_counts.get(src_label, 0) + 1
                    write_manifest_row(mout, ManifestRow(
                        seqid=internal_id,
                        source=src_label,
                        protein_id=prot_id,
                        fetch_type="NCBI_proteome",
                        source_accession=prot_id,
                        returned_fasta_id=prot_id,
                        assembly_or_uniprot=asm,
                        ncbi_taxid=taxid,
                        genus=genus,
                        species=species_only,
                        strain_or_gene=gene,
                        gene_name=gene,
                        protein_name=protein_name,
                        description=protein_name,
                    ))

                    refseq_written += 1

            if manifest_rows and refseq_written == 0:
                # OK rows existed but none yielded a protein -> genuine problem.
                raise RuntimeError("HARD FAIL: wrote 0 proteins from proteomes_manifest.")

            if not manifest_rows:
                # Header-only / all-failed manifest: not an error (see
                # load_proteomes_manifest). Build continues from the other
                # (mod_II / mod_IV / user) sources.
                info("No usable proteome-manifest rows; no proteome-derived proteins added.")
            else:
                info(f"Manifest-derived proteins written: {refseq_written}")
            source_counts["_manifest_total"] = refseq_written

        # ---------------------------------------------------
        # USER-SUPPLIED PROTEINS
        # ---------------------------------------------------
        if args.user_fasta is not None:
            user_meta = {}

            if args.user_metadata is not None:
                user_meta = load_user_metadata(args.user_metadata)

            # Auto-build metadata for fetched proteins from the FASTA itself, so
            # a bare accession list needs no metadata file. Fallback-only:
            # anything supplied above is preserved.
            seed_fetched_protein_stubs(user_meta, args.user_fasta)

            if user_meta and args.ncbi_email:
                user_meta = enrich_user_metadata(
                    user_meta,
                    ncbi_email=args.ncbi_email,
                    api_key=args.ncbi_api_key,
                )
            elif user_meta and not args.ncbi_email:
                warn(
                    "No --ncbi_email supplied; skipping NCBI taxonomy enrichment. "
                    "Genus/species and taxid used exactly as provided / harvested."
                )

            user_written       = 0
            user_source_counts = {}
            user_meta_unused   = set(user_meta.keys())

            for header, seq in iter_fasta(args.user_fasta):
                parts      = header.split()
                pid        = parts[0] if parts else "NA"
                fasta_desc = " ".join(parts[1:]) if len(parts) > 1 else ""
                header_kv  = parse_fasta_header_keyvals(header)

                if pid == "NA":
                    raise RuntimeError(
                        f"HARD FAIL: could not parse protein_id from user FASTA "
                        f"(header='{header[:80]}')"
                    )

                if user_meta:
                    if pid not in user_meta:
                        raise RuntimeError(
                            f"HARD FAIL: no exact metadata row for user protein '{pid}'. "
                            "maniFasta does not infer aliases or guess identifier equivalence. "
                            "Make the FASTA first token exactly match metadata protein_id."
                        )
                    m = user_meta[pid]
                    user_meta_unused.discard(pid)
                else:
                    m = dict(USER_DEFAULTS)

                source         = m["source"]
                taxid          = m["ncbi_taxid"]
                genus          = m["genus"]
                species        = m["species"]
                strain_or_gene = m["strain_or_gene"]
                assembly       = m["assembly"]

                fetch_id = choose_provenance_value(
                    m.get("fetch_id"), header_kv.get("fetch_id"), default="NA"
                )
                fetch_type = choose_provenance_value(
                    m.get("fetch_type"), header_kv.get("fetch_type"), default="user_fasta"
                )
                source_accession = choose_provenance_value(
                    m.get("source_accession"), header_kv.get("source_accession"), default=pid
                )
                returned_fasta_id = choose_provenance_value(
                    m.get("returned_fasta_id"), header_kv.get("returned_fasta_id"), default=pid
                )

                description = (
                    m["description"]
                    if m["description"] != "NA"
                    else (fasta_desc or "NA")
                )

                organism = f"{genus} {species}".strip()

                if taxid != USER_DEFAULTS["ncbi_taxid"] and (genus == "NA" or species == "NA"):
                    raise RuntimeError(
                        f"HARD FAIL: taxid supplied but genus or species is NA for "
                        f"user protein '{pid}' even after NCBI enrichment. "
                        f"Please supply both genus and species in the metadata table."
                    )

                internal_id = next_internal_id(DB_PREFIX, uid)
                uid += 1

                # Module 3 gene/protein-name resolution, NA when unsupplied:
                #   gene         : explicit metadata 'gene' > FASTA-header GN=
                #                  token > legacy strain_or_gene > NA
                #   protein_name : explicit metadata 'protein_name'
                #                  > metadata/FASTA description (cleaned) > NA
                gene = choose_provenance_value(
                    m.get("gene"), header_kv.get("GN"), strain_or_gene, default="NA"
                )
                protein_name = choose_provenance_value(
                    m.get("protein_name"),
                    clean_protein_name(description),
                    default="NA",
                )

                out_header = format_defline(
                    internal_id, source, pid, organism, taxid,
                    assembly=assembly, asm_src=fetch_type,   # fetch_type doubles as src=
                    gene=gene, protein_name=protein_name,
                )

                write_fasta_record(fout, out_header, seq)

                write_manifest_row(mout, ManifestRow(
                    seqid=internal_id,
                    source=source,
                    protein_id=pid,
                    fetch_id=fetch_id,
                    fetch_type=fetch_type,
                    source_accession=source_accession,
                    returned_fasta_id=returned_fasta_id,
                    assembly_or_uniprot=assembly,
                    ncbi_taxid=taxid,
                    genus=genus,
                    species=species,
                    strain_or_gene=strain_or_gene,
                    gene_name=gene,
                    protein_name=protein_name,
                    description=description,
                ))

                user_written += 1
                source_counts[source] = source_counts.get(source, 0) + 1
                user_source_counts[source] = user_source_counts.get(source, 0) + 1

            if user_written == 0:
                raise RuntimeError(
                    f"HARD FAIL: wrote 0 proteins from user FASTA: {args.user_fasta}"
                )

            if user_meta_unused:
                # Metadata rows for accessions a fetcher skipped (--allow-missing)
                # legitimately have no FASTA record; excuse them rather than abort.
                truly_unused = user_meta_unused - skipped_ids
                if truly_unused:
                    raise RuntimeError(
                        f"HARD FAIL: {len(truly_unused)} metadata rows had no exact matching FASTA "
                        f"entry. First unmatched metadata keys: {sorted(truly_unused)[:50]}"
                    )
                warn(
                    f"{len(user_meta_unused)} metadata rows matched skipped/missing accessions; "
                    f"not treated as an error under --allow-missing."
                )

            info(f"User proteins written: {user_written}; exact metadata matches: {user_written if user_meta else 0}")
            if user_source_counts:
                info(
                    "User protein source breakdown: "
                    + ", ".join(
                        f"{label}={count}"
                        for label, count in sorted(user_source_counts.items())
                    )
                )
            source_counts["_user_total"] = user_written

    total = uid - 1
    info(f"Done. Final number of records written: {total}")
    info(f"Output FASTA: {args.out_fasta}")
    info(f"Output mapping: {args.out_map}")

    if args.out_summary is not None:
        write_build_summary(
            Path(args.out_summary),
            args           = args,
            source_counts  = source_counts,
            total          = total,
            start_time     = start_time,
            run_label      = getattr(args, "run_label", None) or "run",
            failed_sources = failed_sources,
        )


# -----------------------------------------------------------
# ARGUMENT PARSING
# -----------------------------------------------------------

def main():
    global DB_PREFIX
    parser = argparse.ArgumentParser(
        description=(
            "Build combined FASTA and mapping file using "
            "HOMD (optional, filterable by body site) + human UniProt-like FASTA "
            "+ SARS2 + cRAP + downloaded proteomes + optional user-supplied proteins."
        )
    )

    # HOMD — all three are optional but must be supplied together
    homd_grp = parser.add_argument_group(
        "HOMD (optional)",
        "All three HOMD arguments must be supplied together, or all omitted. "
        "Pass a body-site-filtered taxonomy TSV (from homd_filter.py) as "
        "--homd_tax to restrict which taxa are included."
    )
    homd_grp.add_argument(
        "--homd_fasta", required=False, type=Path, default=None,
        help="HOMD combined protein FASTA.",
    )
    homd_grp.add_argument(
        "--gca_info", required=False, type=Path, default=None,
        help="GCA_ID_info.txt mapping assembly IDs to HMT IDs.",
    )
    homd_grp.add_argument(
        "--homd_tax", required=False, type=Path, default=None,
        help=(
            "HOMD taxonomy TSV (full or body-site-filtered). "
            "HMT IDs absent from this file are silently excluded from the build."
        ),
    )

    parser.add_argument("--human_fasta", required=False, type=Path, default=None)

    parser.add_argument(
        "--crap_fasta", required=False, type=Path,
        help="Optional cRAP/contaminants FASTA. Records missing OS/OX are skipped.",
    )

    parser.add_argument(
        "--proteomes_manifest", required=False, type=Path,
        help="proteomes_manifest.tsv. Supports old and new manifest column names.",
    )

    parser.add_argument(
        "--user_fasta", required=False, type=Path,
        help="Optional user-supplied protein FASTA file.",
    )
    parser.add_argument(
        "--user_metadata", required=False, type=Path,
        help=(
            "Optional TSV metadata table for user-supplied proteins. "
            "Required column: protein_id. "
            "Optional columns: source, ncbi_taxid, genus, species, "
            "strain_or_gene, assembly, description."
        ),
    )

    parser.add_argument(
        "--failed_accessions", required=False, type=Path, default=None,
        help=(
            "Optional TSV (source_label, accession, fetch_type, reason) of accessions a "
            "fetcher skipped under --allow-missing. These are surfaced in the build summary "
            "and excused from the user-metadata unused-rows check."
        ),
    )

    parser.add_argument(
        "--ncbi_email", required=False, type=str, default=None,
        help=(
            "E-mail address for NCBI eUtils (required by NCBI policy). "
            "Enables automatic taxid <-> genus/species enrichment of user metadata."
        ),
    )
    parser.add_argument(
        "--ncbi_api_key", required=False, type=str, default=None,
        help="Optional NCBI API key (10 req/sec vs 3 req/sec without key).",
    )
    parser.add_argument("--out_fasta",    required=True,  type=Path)
    parser.add_argument("--out_map",     required=True,  type=Path)
    parser.add_argument("--out_summary", required=False, type=str, default=None,
        help="Optional: write build summary TSV to this path.")
    parser.add_argument("--run_label",   required=False, type=str, default=None,
        help="Run label (e.g. AllOralsDB_v1.0) for the summary file.")
    parser.add_argument("--run_dir",     required=False, type=str, default=None,
        help="Run output directory path for the summary file.")
    parser.add_argument("--db_prefix",   required=False, type=str, default=DB_PREFIX,
        help=("Prefix applied to every internal protein ID / DB header "
              "(e.g. %(default)s-00000000001). Set via db_prefix in "
              "000.maniFasta.config; defaults to %(default)s if omitted."))
    parser.add_argument(
        "--source_plan", required=False, type=Path, default=None,
        help=(
            "Optional: source_plan_normalized.tsv. Enables per-source "
            "lineage_taxid_override=<taxid> options (declared in the "
            "registry's options column) to populate the manifest's "
            "lineage_taxid column. Without this, lineage_taxid simply "
            "mirrors each row's true ncbi_taxid."
        ),
    )

    args = parser.parse_args()

    # Allow the DB / header prefix to be set from the config file
    # (db_prefix -> --db_prefix). Fall back to the module default when the
    # flag is empty or unset so internal IDs never start with a bare "-".
    if args.db_prefix and args.db_prefix.strip():
        DB_PREFIX = args.db_prefix.strip()

    build_db(args)


if __name__ == "__main__":
    main()
