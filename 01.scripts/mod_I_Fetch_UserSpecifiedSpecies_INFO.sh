#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# mod_I_Fetch_UserSpecifiedSpecies_INFO.sh  v3.0
#
# Resolve species from a unified input TSV.
# Accepts any combination of columns:
#   Species name    Accession    TaxonID
#
# Priority per row (first match wins):
#   PATH A — Accession present
#     Fetch GenBank record → confirm taxid/organism
#     If Species name absent, use organism from GenBank record
#     download_plan = OK_ACCESSION_PROTEINS
#
#   PATH B — No accession, TaxonID present
#     Look up scientific name from taxid via NCBI efetch
#     download_plan = OK_DATASETS
#
#   PATH C — No accession, no taxid, Species name present
#     NCBI Taxonomy esearch by name → taxid
#     download_plan = OK_DATASETS
#
#   Rows with all three fields empty are skipped.
#   At least one of 'Species name', 'Accession', or 'TaxonID'
#   must be present as a column.
#
# Output columns:
#   Species_name | taxid | organism | input_accession |
#   resolve_status | taxid_source | download_plan
# ============================================================

usage() {
  cat <<'USAGE'
Usage:
  bash mod_I_Fetch_UserSpecifiedSpecies_INFO.sh \
    -i  Species_fetch_list.tsv  \
    -o  species_resolved.tsv    \
    [--sleep 2.0]               \
    [--email you@uni.edu]

Input TSV must have at least one of:
  'Accession'    - GenBank/RefSeq accession  (PATH A, highest priority)
  'TaxonID'      - NCBI taxonomy ID          (PATH B)
  'Species name' - scientific name           (PATH C, lowest priority)

Rows may have any combination of columns.
Priority: Accession > TaxonID > Species name
USAGE
}

INPUT=""
OUTPUT=""
SLEEP_SEC="2.0"
EMAIL=""
API_KEY=""
SOURCE_LABEL=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   INPUT="$2";     shift 2 ;;
    -o|--output)  OUTPUT="$2";    shift 2 ;;
    --sleep)      SLEEP_SEC="$2"; shift 2 ;;
    --email)      EMAIL="$2";     shift 2 ;;
    --api_key)      API_KEY="$2";      shift 2 ;;
    --source_label) SOURCE_LABEL="$2"; shift 2 ;;
    -h|--help)      usage; exit 0 ;;
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
TOOLNAME="maniFasta_step1_v3"
COMMON="tool=${TOOLNAME}"

if [[ -n "$EMAIL" ]]; then
  COMMON="${COMMON}&email=$(python3 -c \
    "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))" "$EMAIL")"
fi
if [[ -n "$API_KEY" ]]; then
  COMMON="${COMMON}&api_key=${API_KEY}"
fi

# ── Helpers ──────────────────────────────────────────────────

sleep_between() {
  python3 -c "import time; time.sleep(float('$SLEEP_SEC'))"
}

fetch_to_file() {
  local url="$1" out="$2"
  for attempt in 1 2 3 4 5; do
    if curl --http1.1 -g -sS -L \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 20 --max-time 120 \
        "$url" -o "$out" && [[ -s "$out" ]]; then
      return 0
    fi
    echo "WARN: fetch failed (attempt $attempt): $url" >&2
    sleep_between
  done
  return 1
}

# Only succeeds if response starts with LOCUS (valid GenBank flat file)
fetch_gb_to_file() {
  local id="$1" out="$2"
  local url="${EUTILS_BASE}/efetch.fcgi?db=nuccore&id=${id}&rettype=gb&retmode=text&${COMMON}"
  for attempt in 1 2 3 4 5; do
    if curl --http1.1 -g -sS -L \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 20 --max-time 120 \
        "$url" -o "$out" \
        && [[ -s "$out" ]] \
        && head -n1 "$out" | grep -q '^LOCUS'; then
      return 0
    fi
    echo "WARN: GB fetch failed/invalid (attempt $attempt): $id" >&2
    sleep_between
  done
  return 1
}

parse_taxid_from_gb()    { grep -m1 -o 'taxon:[0-9]\+' "$1" | cut -d: -f2 || true; }
parse_organism_from_gb() {
  awk '/^  ORGANISM/{sub(/^  ORGANISM[[:space:]]+/,""); print; exit}' "$1" || true
}

get_taxid_from_taxonomy() {
  local spp="$1"
  local term tmp uid

  term="$(python3 - "$spp" <<'PY'
import urllib.parse, sys
s = sys.argv[1].strip()
q = f'"{s}"[Scientific Name] OR "{s}"[All Names]'
print(urllib.parse.quote(q))
PY
)"

  tmp="$(mktemp -p "$TMPDIR" tax_esearch.XXXXXX.xml)"
  if ! fetch_to_file \
      "${EUTILS_BASE}/esearch.fcgi?db=taxonomy&term=${term}&retmode=xml&retmax=1&${COMMON}" \
      "$tmp"; then
    rm -f "$tmp"; return 1
  fi

  uid="$(sed -n 's/.*<Id>\([0-9]\+\)<\/Id>.*/\1/p' "$tmp" | head -n1 || true)"
  rm -f "$tmp"
  [[ -n "$uid" ]] || return 1
  printf '%s\n' "$uid"
}

get_taxonomy_scientific_name() {
  local taxid="$1"
  local tmp name attempt

  # A curl/fetch_to_file "success" (HTTP 200, non-empty body) does not
  # guarantee a real TaxaSet response -- NCBI can return a transient
  # rate-limit/error body that still satisfies that check. A genuinely
  # invalid taxid still returns a body, just without a <ScientificName>, so
  # we cannot tell the two apart from one attempt. Retry with backoff before
  # treating an empty name as confirmed, rather than failing permanently on
  # the first miss.
  for attempt in 1 2 3; do
    tmp="$(mktemp -p "$TMPDIR" tax_efetch.XXXXXX.xml)"
    if fetch_to_file \
        "${EUTILS_BASE}/efetch.fcgi?db=taxonomy&id=${taxid}&retmode=xml&${COMMON}" \
        "$tmp"; then
      name="$(sed -n 's:.*<ScientificName>\(.*\)</ScientificName>.*:\1:p' "$tmp" | head -n1 || true)"
      rm -f "$tmp"
      if [[ -n "$name" ]]; then
        printf '%s\n' "$name"
        return 0
      fi
      echo "   WARN: efetch returned no ScientificName for taxid=${taxid} (attempt ${attempt}/3); retrying" >&2
    else
      rm -f "$tmp"
    fi
    sleep "$((attempt * 2))"
  done
  return 1
}

# ── Assembly metadata from ncbi-datasets ─────────────────────────────────────
# Called when PATH A detects a GCF/GCA accession.
# Uses 'datasets summary genome accession' to get taxid, organism, strain.
# Returns tab-separated: taxid\torganism\tstrain

get_assembly_metadata() {
  local acc="$1"
  local json taxid organism strain

  json="$(datasets summary genome accession "$acc" 2>/dev/null || true)"

  if [[ -z "$json" ]]; then
    echo "" ; return 1
  fi

  taxid="$(printf '%s' "$json" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    rpts = d.get('reports', [])
    if rpts:
        print(rpts[0].get('organism',{}).get('tax_id',''))
except: pass
" 2>/dev/null || true)"

  organism="$(printf '%s' "$json" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    rpts = d.get('reports', [])
    if rpts:
        print(rpts[0].get('organism',{}).get('organism_name',''))
except: pass
" 2>/dev/null || true)"

  strain="$(printf '%s' "$json" | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    rpts = d.get('reports', [])
    if rpts:
        inf = rpts[0].get('assembly_info',{})
        print(inf.get('biosample',{}).get('attributes',[{}])[0].get('value','') or
              rpts[0].get('organism',{}).get('infraspecific_names',{}).get('strain',''))
except: pass
" 2>/dev/null || true)"

  printf '%s\t%s\t%s\n' "${taxid:-}" "${organism:-}" "${strain:-}"
}

# ── PATH A fallback: if accession fails, try taxid then species name ──────────

path_a_fallback() {
  local spp="$1" taxid_hint="$2" acc="$3" fail_reason="$4"
  echo "   WARN: ${fail_reason} for acc=${acc}" >&2

  # Try taxid hint first (PATH B logic)
  if [[ -n "$taxid_hint" ]]; then
    echo "   Falling back to supplied TaxonID: ${taxid_hint}" >&2
    local org
    org="$(get_taxonomy_scientific_name "$taxid_hint" || true)"
    sleep_between
    [[ -n "$org" ]] || org="${spp:-UNKNOWN}"
    local label="${spp:-$org}"
    echo "   taxid=$taxid_hint  organism=$org  (taxid fallback)" >&2
    echo -e "${source_label_row:-UNKNOWN}\t${label}\t${taxid_hint}\t${org}\t${acc}\tOK_FALLBACK\tinput_taxid\tOK_DATASETS" \
      >> "$OUTPUT"
    return
  fi

  # Then try species name (PATH C logic)
  if [[ -n "$spp" ]]; then
    echo "   Falling back to taxonomy search for: ${spp}" >&2
    local taxid org
    taxid="$(get_taxid_from_taxonomy "$spp" || true)"
    sleep_between

    if [[ -z "$taxid" ]]; then
      echo "   WARN: all fallbacks failed for: ${spp}" >&2
      echo -e "${source_label_row:-UNKNOWN}\t${spp}\t\t\t${acc}\tFAILED_ALL\ttaxonomy_fallback\tFAILED_RESOLVE" \
        >> "$OUTPUT"
      return
    fi

    org="$(get_taxonomy_scientific_name "$taxid" || true)"
    sleep_between
    [[ -n "$org" ]] || org="$spp"

    echo "   taxid=$taxid  organism=$org  (name fallback)" >&2
    echo -e "${source_label_row:-UNKNOWN}\t${spp}\t${taxid}\t${org}\t${acc}\tOK_FALLBACK\ttaxonomy_fallback\tOK_DATASETS" \
      >> "$OUTPUT"
    return
  fi

  # Nothing left to try
  echo "   WARN: no fallback available (no taxid or species name supplied)" >&2
  echo -e "${source_label_row:-UNKNOWN}\t\t\t\t${acc}\tFAILED_ALL\tnone\tFAILED_RESOLVE" >> "$OUTPUT"
}

# ── Parse header indices ──────────────────────────────────────
# Accession, TaxonID, and Species name are all optional columns.
# At least one must be present; rows with all fields empty are skipped.

read -r IDX_SPP IDX_ACC IDX_TAX IDX_LABEL < <(python3 - "$INPUT" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8") as f:
    h = next(csv.reader(f, delimiter="\t"))
# Strip comment lines before checking header
h = [c.strip() for c in h]
# Accept both old and new column name variants
SPP_NAMES   = ["species", "Species name", "species name"]
ACC_NAMES   = ["accession", "Accession"]
TAX_NAMES   = ["taxonID", "TaxonID", "taxid"]
LABEL_NAMES = ["source_label"]

h_lower = [x.lower().strip() for x in h]
idx_spp   = next((h_lower.index(n.lower()) for n in SPP_NAMES   if n.lower() in h_lower), -1)
idx_acc   = next((h_lower.index(n.lower()) for n in ACC_NAMES   if n.lower() in h_lower), -1)
idx_tax   = next((h_lower.index(n.lower()) for n in TAX_NAMES   if n.lower() in h_lower), -1)
idx_label = next((h_lower.index(n.lower()) for n in LABEL_NAMES if n.lower() in h_lower), -1)

if idx_spp == -1 and idx_acc == -1 and idx_tax == -1:
    print(
        "ERROR: input TSV must have at least one of: "
        "'species', 'accession', 'taxonID' (case-insensitive). "
        f"Found columns: {h}",
        file=sys.stderr
    )
    sys.exit(2)
print(idx_spp, idx_acc, idx_tax, idx_label)
PY
)

echo "Column indices -> 'species': $IDX_SPP  |  'accession': $IDX_ACC  |  'taxonID': $IDX_TAX  |  'source_label': $IDX_LABEL  (-1 = absent)" >&2

# ── Output header ─────────────────────────────────────────────

echo -e "source_label\tSpecies_name\ttaxid\torganism\tinput_accession\tresolve_status\ttaxid_source\tdownload_plan" \
  > "$OUTPUT"

# ── Main loop ─────────────────────────────────────────────────
# awk -F'\t' splits on literal tab, preserving empty fields correctly.
# It emits: spp|acc|taxid on each data line using | as a safe delimiter,
# skipping the header (NR==1) and comment lines.
# The while loop then splits on | — a non-whitespace char, no IFS collapsing.

while IFS='|' read -r spp acc taxid_in source_label_row; do

  # Skip rows with all fields empty
  [[ -n "$spp" || -n "$acc" || -n "$taxid_in" ]] || continue

  echo "==> spp='${spp:-<empty>}'  acc='${acc:-<none>}'  taxid='${taxid_in:-<none>}'" >&2

  # ════════════════════════════════════════════════════════════
  # PATH A: Accession present (highest priority)
  # ════════════════════════════════════════════════════════════
  if [[ -n "$acc" ]]; then

    # GCF/GCA are assembly accessions — fetch metadata via ncbi-datasets,
    # then write download_plan=OK_ASSEMBLY_ACCESSION so the download step
    # uses that specific assembly rather than auto-picking by taxid.
    if [[ "$acc" =~ ^GC[FA]_ ]]; then
      echo "   [PATH A] assembly accession (${acc}); fetching metadata via datasets" >&2

      meta="$(get_assembly_metadata "$acc" || true)"
      asm_taxid="$(printf '%s' "$meta" | cut -f1)"
      asm_org="$(printf '%s'   "$meta" | cut -f2)"
      asm_strain="$(printf '%s' "$meta" | cut -f3)"

      # -- Assembly accession did not resolve ------------------------------
      # A withdrawn, suppressed, or typo'd GCA/GCF returns nothing from
      # `datasets summary genome accession`. Previously the supplied taxid was
      # substituted in here and the row was STILL written as
      #   resolve_status=OK | taxid_source=input_accession | OK_ASSEMBLY_ACCESSION
      # which (a) recorded a failed lookup as a clean accession-based resolve,
      # (b) mis-attributed the taxid's provenance, and (c) pinned the download
      # step to the dead accession -- it retried that accession 3x instead of
      # reaching the refseq-reference > representative > any > genbank chain.
      #
      # Hand off to the documented accession -> taxid -> species chain instead.
      # path_a_fallback keeps ${acc} in the input_accession column for audit and
      # marks resolve_status=OK_FALLBACK with the true taxid_source, so the
      # downgrade stays greppable in the resolved-species TSV that ships in
      # run_setup/. NOTE: the taxid substitution that used to sit here is
      # deliberately gone -- asm_taxid must reflect the accession lookup ONLY.
      if [[ -z "$asm_taxid" ]]; then
        path_a_fallback "$spp" "$taxid_in" "$acc" "assembly metadata lookup failed"
        sleep_between; continue
      fi

      [[ -n "$asm_org"    ]] || asm_org="${spp}"
      [[ -n "$asm_strain" ]] || asm_strain=""

      # Label: prefer user-supplied species name, fall back to organism from NCBI
      label="${spp:-$asm_org}"

      echo "   taxid=${asm_taxid}  organism=${asm_org}  strain=${asm_strain}  label=${label}" >&2
      echo -e "${source_label_row:-UNKNOWN}\t${label}\t${asm_taxid}\t${asm_org}\t${acc}\tOK\tinput_accession\tOK_ASSEMBLY_ACCESSION" \
        >> "$OUTPUT"
      sleep_between; continue
    fi

    echo "   [PATH A] accession-based fetch" >&2
    tmpgb="$(mktemp -p "$TMPDIR" gb.XXXXXX)"

    if ! fetch_gb_to_file "$acc" "$tmpgb"; then
      rm -f "$tmpgb"
      path_a_fallback "$spp" "$taxid_in" "$acc" "accession fetch failed"
      sleep_between; continue
    fi

    taxid="$(parse_taxid_from_gb  "$tmpgb")"
    org="$(parse_organism_from_gb "$tmpgb")"
    rm -f "$tmpgb"

    # If species name not supplied, use organism from GenBank
    [[ -n "$spp" ]] || spp="$org"

    if [[ -z "$taxid" ]]; then
      path_a_fallback "$spp" "$taxid_in" "$acc" "no taxid in GenBank record"
      sleep_between; continue
    fi

    echo "   taxid=$taxid  organism=$org  label=$spp" >&2
    echo -e "${source_label_row:-UNKNOWN}\t${spp}\t${taxid}\t${org}\t${acc}\tOK\tinput_accession\tOK_ACCESSION_PROTEINS" \
      >> "$OUTPUT"

  # ════════════════════════════════════════════════════════════
  # PATH B: No accession, TaxonID present
  # ════════════════════════════════════════════════════════════
  elif [[ -n "$taxid_in" ]]; then

    echo "   [PATH B] taxid-based lookup" >&2
    org="$(get_taxonomy_scientific_name "$taxid_in" || true)"
    sleep_between

    if [[ -z "$org" ]]; then
      echo "   WARN: could not resolve organism name for taxid=${taxid_in}" >&2
      # If we also have a species name, use it as the label
      local_label="${spp:-UNKNOWN}"
      echo -e "${source_label_row:-UNKNOWN}\t${local_label}\t${taxid_in}\t\t\tFAILED_TAXID_LOOKUP\tinput_taxid\tFAILED_RESOLVE" \
        >> "$OUTPUT"
      sleep_between; continue
    fi

    # Use supplied species name as label if present, else use NCBI name
    label="${spp:-$org}"
    echo "   taxid=$taxid_in  organism=$org  label=$label" >&2
    echo -e "${source_label_row:-UNKNOWN}\t${label}\t${taxid_in}\t${org}\t\tOK\tinput_taxid\tOK_DATASETS" \
      >> "$OUTPUT"

  # ════════════════════════════════════════════════════════════
  # PATH C: No accession, no taxid — species name only (lowest priority)
  # ════════════════════════════════════════════════════════════
  else

    echo "   [PATH C] name-based taxonomy search" >&2
    taxid="$(get_taxid_from_taxonomy "$spp" || true)"
    sleep_between

    if [[ -z "$taxid" ]]; then
      echo "   WARN: no taxid found for: $spp" >&2
      echo -e "${source_label_row:-UNKNOWN}\t${spp}\t\t\t\tFAILED_TAXONOMY_RESOLVE\ttaxonomy_esearch\tFAILED_RESOLVE" \
        >> "$OUTPUT"
      sleep_between; continue
    fi

    org="$(get_taxonomy_scientific_name "$taxid" || true)"
    sleep_between
    [[ -n "$org" ]] || org="$spp"

    echo "   taxid=$taxid  organism=$org" >&2
    echo -e "${source_label_row:-UNKNOWN}\t${spp}\t${taxid}\t${org}\t\tOK\ttaxonomy_esearch\tOK_DATASETS" \
      >> "$OUTPUT"
  fi

  sleep_between
done < <(
  awk -F'\t' -v ispp="$IDX_SPP" -v iacc="$IDX_ACC" -v itax="$IDX_TAX" -v ilabel="$IDX_LABEL" '
    NR == 1 { next }          # skip header
    /^#/    { next }          # skip comment lines
    /^[[:space:]]*$/ { next } # skip blank lines
    {
      spp   = (ispp   >= 0) ? $(ispp+1)   : ""
      acc   = (iacc   >= 0) ? $(iacc+1)   : ""
      tax   = (itax   >= 0) ? $(itax+1)   : ""
      lbl   = (ilabel >= 0) ? $(ilabel+1) : ""
      # skip all-empty rows
      if (spp == "" && acc == "" && tax == "") next
      # trim leading/trailing whitespace
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", spp)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", acc)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", tax)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", lbl)
      print spp "|" acc "|" tax "|" lbl
    }
  ' "$INPUT"
)

echo "Done. Wrote: $OUTPUT" >&2
