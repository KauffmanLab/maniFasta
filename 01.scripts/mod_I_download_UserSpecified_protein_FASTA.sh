#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# mod_I_download_UserSpecified_protein_FASTA.sh
#
# Reads species_resolved.tsv from step 1 and downloads
# protein FASTAs using two strategies:
#
#   download_plan=OK_ACCESSION_PROTEINS
#     -> efetch on input_accession.
#        Accession type is auto-detected:
#          Protein accession  -> db=protein  rettype=fasta
#          Nucleotide accession -> db=nuccore rettype=fasta_cds_aa
#
#   download_plan=OK_DATASETS
#     -> ncbi-datasets-cli, best RefSeq/GenBank assembly by taxid
#        priority: RefSeq reference > representative > any > GenBank
#
# Output per species:
#   <safe_name>__taxid<TAXID>__<acc_or_asm>.protein.faa
# Plus:
#   proteomes_manifest.tsv
# ============================================================

usage() {
  cat <<'USAGE'
Usage:
  ./mod_I_download_UserSpecified_protein_FASTA.sh \
    -i  species_resolved.tsv  \
    -o  /path/to/outdir       \
    [--sleep 1.0]             \
    [--email you@uni.edu]

Requires: curl, python3, datasets (ncbi-datasets-cli), unzip
USAGE
}

INPUT=""
OUTDIR=""
SLEEP_SEC="1.0"
EMAIL=""
API_KEY=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)  INPUT="$2";     shift 2 ;;
    -o|--outdir) OUTDIR="$2";    shift 2 ;;
    --sleep)     SLEEP_SEC="$2"; shift 2 ;;
    --email)     EMAIL="$2";     shift 2 ;;
    --api_key)   API_KEY="$2";   shift 2 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$INPUT"  && -n "$OUTDIR" ]] || { echo "ERROR: -i and -o required" >&2; exit 1; }
[[ -f "$INPUT" ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }

command -v curl     >/dev/null || { echo "ERROR: curl not found"     >&2; exit 1; }
command -v python3  >/dev/null || { echo "ERROR: python3 not found"  >&2; exit 1; }
command -v datasets >/dev/null || { echo "ERROR: datasets (ncbi-datasets-cli) not found" >&2; exit 1; }
command -v unzip    >/dev/null || { echo "ERROR: unzip not found"    >&2; exit 1; }

: "${TMPDIR:=/tmp}"
mkdir -p "$TMPDIR" "$OUTDIR"

EUTILS_BASE="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOLNAME="maniFasta_step2_v2"
COMMON="tool=${TOOLNAME}"

if [[ -n "$EMAIL" ]]; then
  COMMON="${COMMON}&email=$(python3 -c \
    "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1]))" "$EMAIL")"
fi
if [[ -n "$API_KEY" ]]; then
  COMMON="${COMMON}&api_key=${API_KEY}"
fi

# ── Helpers ─────────────────────────────────────────────────

sleep_between() {
  python3 -c "import time; time.sleep(float('$SLEEP_SEC'))"
}

sanitize_name() {
  python3 - "$1" <<'PY'
import re, sys
s = sys.argv[1].strip() if len(sys.argv) > 1 else ""
s = re.sub(r"\s+", "_", s)
s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
print(s)
PY
}

# ── PATH A helper: classify accession type ──────────────────
#
# NCBI accession formats:
#   RefSeq protein   : NP_, WP_, XP_, YP_, AP_, ZP_ + digits
#   GenBank protein  : 3 uppercase letters + digits (e.g. AAO03598)
#   RefSeq nucleotide: NC_, NM_, NR_, NT_, NW_, NG_, NZ_ + digits
#   GenBank nucleotide: 1-2 uppercase letters + digits (e.g. X73487, AF513913)
#
# Returns "protein" or "nuccore".

classify_accession() {
  local acc="$1"
  # RefSeq protein prefixes
  if echo "$acc" | grep -qE '^(NP|WP|XP|YP|AP|ZP)_[0-9]'; then
    echo "protein"; return
  fi
  # GenBank protein: exactly 3 uppercase letters followed by digits
  if echo "$acc" | grep -qE '^[A-Z]{3}[0-9]'; then
    echo "protein"; return
  fi
  # Everything else treated as nucleotide
  echo "nuccore"
}

# ── PATH A helper: fetch proteins by accession ──────────────
#
# Routes to the correct NCBI database and retrieval type:
#   Protein accession    -> db=protein  rettype=fasta
#   Nucleotide accession -> db=nuccore  rettype=fasta_cds_aa

download_proteins_by_accession() {
  local acc="$1" outfile="$2"
  local db rettype url

  local acc_type
  acc_type="$(classify_accession "$acc")"

  if [[ "$acc_type" == "protein" ]]; then
    db="protein"
    rettype="fasta"
  else
    db="nuccore"
    rettype="fasta_cds_aa"
  fi

  echo "   acc_type=${acc_type}  db=${db}  rettype=${rettype}" >&2
  url="${EUTILS_BASE}/efetch.fcgi?db=${db}&id=${acc}&rettype=${rettype}&retmode=text&${COMMON}"

  for attempt in 1 2 3 4 5; do
    rm -f "$outfile"
    if curl --http1.1 -g -sS -L \
        --retry 5 --retry-all-errors --retry-delay 2 \
        --connect-timeout 20 --max-time 180 \
        "$url" -o "$outfile"; then

      # Valid: non-empty file starting with a FASTA header
      if [[ -s "$outfile" ]] && grep -q '^>' "$outfile"; then
        return 0
      fi

      # NCBI sometimes returns an XML error message instead
      if [[ -s "$outfile" ]] && head -n3 "$outfile" | grep -qi 'error\|esummary\|<eFetchResult'; then
        echo "WARN: NCBI returned an error page for acc=$acc:" >&2
        head -n5 "$outfile" >&2
        return 1   # Don't retry XML error — it won't change
      fi

      echo "WARN: ${rettype} returned empty/invalid (attempt $attempt): $acc" >&2
    else
      echo "WARN: curl failed (attempt $attempt): $acc" >&2
    fi
    sleep_between
  done
  return 1
}

# ── PATH B helpers: ncbi-datasets-cli ───────────────────────

try_download() {
  local taxid="$1" zipout="$2"; shift 2
  rm -f "$zipout"
  if datasets download genome taxon "$taxid" \
      --include protein --filename "$zipout" "$@" >/dev/null 2>&1; then
    [[ -s "$zipout" ]] && return 0
  fi
  return 1
}

find_assembly_accession() {
  local root="$1" acc=""
  acc="$(find "$root/ncbi_dataset/data" -maxdepth 1 -type d \
      -name 'GCF_*' -printf '%f\n' 2>/dev/null | head -n1 || true)"
  [[ -z "$acc" ]] && \
  acc="$(find "$root/ncbi_dataset/data" -maxdepth 1 -type d \
      -name 'GCA_*' -printf '%f\n' 2>/dev/null | head -n1 || true)"
  printf '%s' "$acc"
}

find_protein_faa() {
  local root="$1" p=""
  p="$(find "$root/ncbi_dataset/data" -maxdepth 2 -type f \
      -name 'protein.faa' 2>/dev/null | head -n1 || true)"
  if [[ -z "$p" ]]; then
    p="$(find "$root/ncbi_dataset/data" -maxdepth 2 -type f \
        -name 'protein.faa.gz' 2>/dev/null | head -n1 || true)"
    if [[ -n "$p" ]]; then
      gunzip -c "$p" > "${p%.gz}"
      p="${p%.gz}"
    fi
  fi
  printf '%s' "$p"
}

# ── Parse header indices ─────────────────────────────────────

read -r IDX_SPP IDX_TAX IDX_ACC IDX_STAT IDX_PLAN IDX_LABEL < <(python3 - "$INPUT" <<'PY'
import csv, sys
with open(sys.argv[1], newline="", encoding="utf-8") as f:
    h = next(csv.reader(f, delimiter="\t"))
need = ["Species_name", "taxid", "input_accession", "resolve_status", "download_plan"]
for n in need:
    if n not in h:
        print(f"ERROR: missing column '{n}'. Found: {h}", file=sys.stderr)
        sys.exit(2)
idx_label = h.index("source_label") if "source_label" in h else -1
print(h.index("Species_name"), h.index("taxid"), h.index("input_accession"),
      h.index("resolve_status"), h.index("download_plan"), idx_label)
PY
)

manifest="$OUTDIR/proteomes_manifest.tsv"
echo -e "source_label\tSpecies_name\ttaxid\tinput_accession\tpicked_asm_or_acc\tsource\tstatus\tout_faa\tn_proteins\tnote" \
  > "$manifest"

echo "== Step 2: protein download =="
echo "Input:  $INPUT"
echo "Outdir: $OUTDIR"
echo

# ── Main loop ────────────────────────────────────────────────

tail -n +2 "$INPUT" | while IFS= read -r _raw; do
  IFS=$'\x01' read -r -a F <<< "${_raw//$'\t'/$'\x01'}"

  spp="${F[$IDX_SPP]:-}";   spp="${spp%$'\r'}";   spp="$(printf '%s' "$spp"  | sed 's/[[:space:]]*$//')"
  src_lbl="${F[$IDX_LABEL]:-}"; src_lbl="${src_lbl%$'\r'}"; src_lbl="$(printf '%s' "$src_lbl" | xargs)"
  taxid="${F[$IDX_TAX]:-}"; taxid="${taxid%$'\r'}"; taxid="$(printf '%s' "$taxid" | tr -d ' ')"
  acc="${F[$IDX_ACC]:-}";   acc="${acc%$'\r'}";   acc="$(printf '%s' "$acc"  | xargs)"
  rstat="${F[$IDX_STAT]:-}"; rstat="${rstat%$'\r'}"
  plan="${F[$IDX_PLAN]:-}";  plan="${plan%$'\r'}"

  # Skip rows that failed resolution in step 1 -- but still record a
  # manifest row (status=$rstat) so this is visible in build_summary.tsv
  # instead of disappearing silently. Every other failure branch in this
  # script writes a row before continuing; this one previously did not.
  if [[ "$rstat" != "OK" && "$rstat" != "OK_FALLBACK" ]]; then
    echo "SKIP (resolve_status=$rstat): $spp" >&2
    echo -e "${src_lbl}\t${spp}\t${taxid}\t\t\t\t${rstat}\t\t0\tfailed at step-1 taxonomy resolution" \
      >> "$manifest"
    continue
  fi

  safe_spp="$(sanitize_name "$spp")"
  echo "==> $spp  plan=$plan  acc=${acc:-<none>}  taxid=${taxid:-<none>}" >&2

  # ── PATH A: accession-based protein download ────────────
  if [[ "$plan" == "OK_ACCESSION_PROTEINS" ]]; then

    if [[ -z "$acc" ]]; then
      echo "  ERROR: plan=OK_ACCESSION_PROTEINS but input_accession is empty" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t\t\t\tFAILED_MISSING_ACCESSION\t\t0\t" >> "$manifest"
      continue
    fi

    out="$OUTDIR/${safe_spp}__taxid${taxid}__acc${acc}.protein.faa"

    if download_proteins_by_accession "$acc" "$out"; then
      nseq="$(grep -c '^>' "$out" 2>/dev/null || true)"
      echo "  -> acc=$acc  proteins=$nseq" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${acc}\tefetch\tOK\t${out}\t${nseq}\t" \
        >> "$manifest"
    else
      echo "  ERROR: download failed for acc=$acc" >&2
      rm -f "$out"
      echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${acc}\tefetch\tFAILED_DOWNLOAD\t\t0\t" \
        >> "$manifest"
    fi

  # ── PATH B: ncbi-datasets best assembly by taxid ────────
  elif [[ "$plan" == "OK_DATASETS" || "$plan" == "OK_FALLBACK" ]]; then

    if [[ -z "$taxid" ]]; then
      echo "  ERROR: plan=${plan} but taxid is empty" >&2
      echo -e "${spp}\t\t\t\t\tFAILED_MISSING_TAXID\t\t0\t" >> "$manifest"
      continue
    fi

    work="$(mktemp -d -p "$TMPDIR" "ds_${safe_spp}.XXXXXX")"
    zip="$work/datasets.zip"
    unpack="$work/unpacked"
    mkdir -p "$unpack"
    src=""; note=""

    # Priority order: RefSeq reference > representative > any RefSeq > GenBank
    if   try_download "$taxid" "$zip" --assembly-source refseq --reference;     then src="refseq_reference";     note="refseq --reference"
    elif try_download "$taxid" "$zip" --assembly-source refseq --representative; then src="refseq_representative"; note="refseq --representative"
    elif try_download "$taxid" "$zip" --assembly-source refseq;                  then src="refseq_any";            note="refseq any"
    elif try_download "$taxid" "$zip" --assembly-source genbank;                 then src="genbank_fallback";      note="genbank fallback"
    else
      echo "  ERROR: datasets download failed for taxid=$taxid" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t\t\t\tFAILED_DOWNLOAD\t\t0\tdatasets failed all strategies" \
        >> "$manifest"
      rm -rf "$work"; sleep_between; continue
    fi

    if ! unzip -q "$zip" -d "$unpack"; then
      echo "  ERROR: unzip failed" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t\t\t${src}\tFAILED_UNZIP\t\t0\tunzip failed" >> "$manifest"
      rm -rf "$work"; sleep_between; continue
    fi

    asm="$(find_assembly_accession "$unpack")"; [[ -n "$asm" ]] || asm="UNKNOWN_ASM"
    prot="$(find_protein_faa "$unpack")"

    if [[ -z "$prot" ]]; then
      echo "  ERROR: protein.faa not found in dataset package" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t\t${asm}\t${src}\tFAILED_NO_PROTEIN\t\t0\tprotein.faa missing" \
        >> "$manifest"
      rm -rf "$work"; sleep_between; continue
    fi

    out="$OUTDIR/${safe_spp}__taxid${taxid}__${asm}.protein.faa"
    cp -f "$prot" "$out"
    nseq="$(grep -c '^>' "$out" 2>/dev/null || true)"
    echo "  -> $asm ($src)  proteins=$nseq" >&2
    echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${asm}\t${src}\tOK\t${out}\t${nseq}\t${note}" >> "$manifest"
    rm -rf "$work"

  # ── PATH C: specific assembly accession download ───────────
  elif [[ "$plan" == "OK_ASSEMBLY_ACCESSION" ]]; then

    if [[ -z "$acc" ]]; then
      echo "  ERROR: plan=OK_ASSEMBLY_ACCESSION but input_accession is empty" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t\t\t\tFAILED_MISSING_ACCESSION\t\t0\t" >> "$manifest"
      continue
    fi

    work="$(mktemp -d -p "$TMPDIR" "asm_${safe_spp}.XXXXXX")"
    zip="$work/datasets.zip"
    unpack="$work/unpacked"
    mkdir -p "$unpack"

    DOWNLOAD_OK=false
    for attempt in 1 2 3; do
      rm -f "$zip"
      if datasets download genome accession "$acc" \
          --include protein --filename "$zip" >/dev/null 2>&1 \
          && [[ -s "$zip" ]]; then
        DOWNLOAD_OK=true
        break
      fi
      echo "  WARN: assembly download attempt ${attempt} failed for ${acc}" >&2
      sleep_between
    done

    if [[ "$DOWNLOAD_OK" != "true" ]]; then
      echo "  ERROR: assembly download failed for acc=$acc" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${acc}\tncbi_datasets\tFAILED_DOWNLOAD\t\t0\t" >> "$manifest"
      rm -rf "$work"; sleep_between; continue
    fi

    if ! unzip -q "$zip" -d "$unpack"; then
      echo "  ERROR: unzip failed for $acc" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${acc}\tncbi_datasets\tFAILED_UNZIP\t\t0\t" >> "$manifest"
      rm -rf "$work"; sleep_between; continue
    fi

    prot="$(find_protein_faa "$unpack")"

    if [[ -z "$prot" ]]; then
      echo "  ERROR: protein.faa not found for $acc" >&2
      echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${acc}\tncbi_datasets\tFAILED_NO_PROTEIN\t\t0\tprotein.faa missing" >> "$manifest"
      rm -rf "$work"; sleep_between; continue
    fi

    out="$OUTDIR/${safe_spp}__taxid${taxid}__${acc}.protein.faa"
    cp -f "$prot" "$out"
    nseq="$(grep -c '^>' "$out" 2>/dev/null || true)"
    echo "  -> ${acc} (specific assembly)  proteins=$nseq" >&2
    echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t${acc}\tncbi_datasets_accession\tOK\t${out}\t${nseq}\tspecific assembly" >> "$manifest"
    rm -rf "$work"

  else
    echo "  SKIP: unrecognized download_plan='$plan'" >&2
    echo -e "${src_lbl}\t${spp}\t${taxid}\t${acc}\t\t\tSKIPPED_UNKNOWN_PLAN\t\t0\tplan=$plan" >> "$manifest"
    continue
  fi

  sleep_between
done

echo
echo "DONE. Manifest: $manifest" >&2
