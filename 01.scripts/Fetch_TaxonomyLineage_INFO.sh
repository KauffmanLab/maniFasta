#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Fetch_TaxonomyLineage_INFO.sh  v1.0
#
# Given a list of unique NCBI taxids, batch-fetches full taxonomic
# lineage via NCBI Taxonomy epost+efetch and writes a lineage cache TSV.
#
# Handles cellular AND viral taxa uniformly:
#   - Cellular lineages carry a 'superkingdom' rank (Bacteria/Archaea/Eukaryota)
#   - Viral lineages carry a 'realm' rank instead (e.g. Riboviria) -- many
#     older/unclassified viral taxa (satellites, unclassified ssRNA, etc.)
#     have NEITHER and will show blank at the top -- this is correct,
#     not a failure; they still roll up at whatever rank IS assigned.
#   - 'superkingdom' and 'realm' are merged into one 'domain_or_realm'
#     column so cellular and viral rows sit in the same column.
#
# Also writes a full root-to-tip name path per taxid (pipe-delimited),
# independent of named ranks -- this is the form Krona and a Sankey
# edge-list both actually want, and it degrades gracefully for sparse
# viral lineages where named ranks run out early.
#
# Output columns:
#   taxid  own_rank  own_name  domain_or_realm  kingdom  phylum  class
#   order  family  genus  species  full_lineage_path
# ============================================================

usage() {
  cat <<'USAGE'
Usage:
  bash Fetch_TaxonomyLineage_INFO.sh \
    -i  unique_taxids.txt        \
    -o  taxid_lineage_cache.tsv  \
    [--batch-size 400]           \
    [--sleep 0.5]                \
    [--email you@uni.edu]        \
    [--api_key XXXX]

Input: one taxid per line (blank lines and non-numeric lines are skipped).
USAGE
}

INPUT=""
OUTPUT=""
BATCH_SIZE="400"
SLEEP_SEC="0.5"
EMAIL=""
API_KEY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)      INPUT="$2";      shift 2 ;;
    -o|--output)     OUTPUT="$2";     shift 2 ;;
    --batch-size)    BATCH_SIZE="$2"; shift 2 ;;
    --sleep)         SLEEP_SEC="$2";  shift 2 ;;
    --email)         EMAIL="$2";      shift 2 ;;
    --api_key)       API_KEY="$2";    shift 2 ;;
    -h|--help)       usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT"  ]] || { echo "ERROR: -i/--input required"  >&2; exit 1; }
[[ -n "$OUTPUT" ]] || { echo "ERROR: -o/--output required" >&2; exit 1; }
[[ -f "$INPUT"  ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

command -v curl    >/dev/null || { echo "ERROR: curl not found"    >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not found" >&2; exit 1; }

: "${TMPDIR:=/tmp}"
mkdir -p "$TMPDIR"

EUTILS_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOLNAME="maniFasta_taxonomyLineage_v1"
COMMON="tool=${TOOLNAME}"

if [[ -n "$EMAIL" ]]; then
  COMMON="${COMMON}&email=$(python3 -c \
    "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))" "$EMAIL")"
fi
if [[ -n "$API_KEY" ]]; then
  COMMON="${COMMON}&api_key=${API_KEY}"
fi

sleep_between() {
  python3 -c "import time; time.sleep(float('$SLEEP_SEC'))"
}

# ── epost: post a batch of taxids, return WebEnv + QueryKey ───────
epost_taxids() {
  local id_csv="$1" out="$2"
  for attempt in 1 2 3 4 5; do
    if curl --http1.1 -g -sS -L \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 20 --max-time 120 \
        --data-urlencode "db=taxonomy" \
        --data-urlencode "id=${id_csv}" \
        "${EUTILS_BASE}/epost.fcgi?${COMMON}" -o "$out" \
        && [[ -s "$out" ]] \
        && grep -q '<WebEnv>' "$out"; then
      return 0
    fi
    echo "WARN: epost failed (attempt $attempt)" >&2
    sleep_between
  done
  return 1
}

# ── efetch: pull full TaxaSet XML using WebEnv/QueryKey ────────────
efetch_lineage_xml() {
  local webenv="$1" querykey="$2" out="$3"
  local url="${EUTILS_BASE}/efetch.fcgi?db=taxonomy&WebEnv=${webenv}&query_key=${querykey}&retmode=xml&${COMMON}"
  for attempt in 1 2 3 4 5; do
    if curl --http1.1 -g -sS -L \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 20 --max-time 180 \
        "$url" -o "$out" \
        && [[ -s "$out" ]] \
        && grep -q '<TaxaSet' "$out"; then
      return 0
    fi
    echo "WARN: efetch failed (attempt $attempt)" >&2
    sleep_between
  done
  return 1
}

# ── parse epost response for WebEnv/QueryKey ───────────────────────
parse_webenv_querykey() {
  local f="$1"
  python3 - "$f" <<'PY'
import re, sys
text = open(sys.argv[1], encoding="utf-8", errors="replace").read()
we = re.search(r"<WebEnv>(.*?)</WebEnv>", text)
qk = re.search(r"<QueryKey>(.*?)</QueryKey>", text)
if not (we and qk):
    sys.exit(1)
print(we.group(1))
print(qk.group(1))
PY
}

# ── parse TaxaSet XML into cache rows ───────────────────────────────
# Merges 'superkingdom' (cellular) and 'realm' (viral) into one
# domain_or_realm column. Also emits the full root-to-tip name path
# (independent of named ranks) for Krona/Sankey use downstream.
parse_taxaset_xml() {
  local f="$1"
  python3 - "$f" <<'PY'
import sys
import xml.etree.ElementTree as ET

RANK_COLS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

try:
    tree = ET.parse(sys.argv[1])
except ET.ParseError as e:
    print(f"PARSE_ERROR: {e}", file=sys.stderr)
    sys.exit(1)

root = tree.getroot()
for taxon in root.findall("Taxon"):
    taxid = taxon.findtext("TaxId", "")
    name = taxon.findtext("ScientificName", "")
    own_rank = taxon.findtext("Rank", "no rank")

    rank_map = {}
    path = []
    lineage_ex = taxon.find("LineageEx")
    if lineage_ex is not None:
        for t in lineage_ex.findall("Taxon"):
            r = t.findtext("Rank", "no rank")
            n = t.findtext("ScientificName", "")
            if n:
                path.append(n)
            if r == "superkingdom" or r == "realm":
                rank_map.setdefault("domain_or_realm", n)
            elif r in RANK_COLS:
                rank_map[r] = n

    # include the queried taxon itself in the path and rank map
    if name:
        path.append(name)
    if own_rank == "superkingdom" or own_rank == "realm":
        rank_map.setdefault("domain_or_realm", name)
    elif own_rank in RANK_COLS:
        rank_map[own_rank] = name

    row = [
        taxid, own_rank, name,
        rank_map.get("domain_or_realm", ""),
    ] + [rank_map.get(r, "") for r in RANK_COLS] + [
        "|".join(path)
    ]
    print("\t".join(row))
PY
}

# ── Build deduped, validated taxid list ─────────────────────────────
clean_list="$(mktemp -p "$TMPDIR" taxids_clean.XXXXXX)"
grep -E '^[0-9]+$' "$INPUT" | sort -u > "$clean_list"
n_total="$(grep -c . "$clean_list" || true)"
n_skipped=$(( $(grep -c . "$INPUT" || true) - n_total ))
echo "Unique numeric taxids: $n_total  (skipped non-numeric/blank: $n_skipped)" >&2

echo -e "taxid\town_rank\town_name\tdomain_or_realm\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tfull_lineage_path" \
  > "$OUTPUT"

# ── Batch loop ────────────────────────────────────────────────────
batch_num=0
n_resolved=0
failed_log="${OUTPUT%.tsv}.failed_batches.log"
: > "$failed_log"

while IFS= read -r -d '' batch_ids; do
  batch_num=$((batch_num + 1))
  n_ids="$(printf '%s' "$batch_ids" | tr ',' '\n' | grep -c .)"
  echo "[batch $batch_num] $n_ids taxids" >&2

  epost_out="$(mktemp -p "$TMPDIR" epost.XXXXXX.xml)"
  if ! epost_taxids "$batch_ids" "$epost_out"; then
    echo "ERROR: epost failed for batch $batch_num" >&2
    echo "epost_failed: $batch_ids" >> "$failed_log"
    rm -f "$epost_out"; continue
  fi

  read -r webenv querykey < <(parse_webenv_querykey "$epost_out" | tr '\n' ' ') || {
    echo "ERROR: could not parse WebEnv/QueryKey for batch $batch_num" >&2
    echo "parse_failed: $batch_ids" >> "$failed_log"
    rm -f "$epost_out"; continue
  }
  rm -f "$epost_out"

  efetch_out="$(mktemp -p "$TMPDIR" efetch.XXXXXX.xml)"
  if ! efetch_lineage_xml "$webenv" "$querykey" "$efetch_out"; then
    echo "ERROR: efetch failed for batch $batch_num" >&2
    echo "efetch_failed: $batch_ids" >> "$failed_log"
    rm -f "$efetch_out"; continue
  fi

  n_before="$(grep -c . "$OUTPUT")"
  parse_taxaset_xml "$efetch_out" >> "$OUTPUT" || {
    echo "ERROR: XML parse failed for batch $batch_num" >&2
    echo "xml_parse_failed: $batch_ids" >> "$failed_log"
  }
  n_after="$(grep -c . "$OUTPUT")"
  n_resolved=$((n_resolved + (n_after - n_before)))
  rm -f "$efetch_out"

  sleep_between
done < <(
  awk -v bs="$BATCH_SIZE" '
    { ids[NR]=$0 }
    END {
      n = NR; buf = ""
      for (i=1; i<=n; i++) {
        buf = (buf=="") ? ids[i] : buf "," ids[i]
        if (i % bs == 0 || i == n) {
          printf "%s\0", buf
          buf = ""
        }
      }
    }
  ' "$clean_list"
)

echo >&2
echo "Done. Requested taxids: $n_total | Lineage rows written: $n_resolved" >&2
echo "Output: $OUTPUT" >&2
if [[ -s "$failed_log" ]]; then
  echo "Some batches failed -- see: $failed_log" >&2
else
  rm -f "$failed_log"
fi
